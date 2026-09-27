"""Loopback-only HTTP server; no additional web framework or frontend build."""
from __future__ import annotations

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib.resources import files
import json
from pathlib import Path
import secrets
from urllib.parse import parse_qs, quote, urlsplit
import webbrowser

from .service import WebService


class WebServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self, service: WebService, port=8765):
        self.service = service
        self.token = secrets.token_urlsafe(32)
        super().__init__(("127.0.0.1", port), Handler)
        self.allowed_hosts = {f"127.0.0.1:{self.server_port}", f"localhost:{self.server_port}"}


class Handler(BaseHTTPRequestHandler):
    def log_message(self, format, *args):
        # Keep polling quiet and never log form contents or request tokens.
        pass

    def send_content(self, body: bytes, content_type="application/json; charset=utf-8", status=200, download=None):
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy", "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'; base-uri 'none'; form-action 'self'")
        if download:
            self.send_header("Content-Disposition", f'attachment; filename="{download}"')
        self.end_headers()
        self.wfile.write(body)

    def send_path(self, path: Path, content_type: str, download: str):
        size = path.stat().st_size
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(size))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Disposition", "attachment; filename*=UTF-8''" + quote(Path(download).name))
        self.end_headers()
        with path.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                self.wfile.write(chunk)

    def send_json(self, value, status=200):
        self.send_content(json.dumps(value, ensure_ascii=False, allow_nan=False).encode("utf-8"), status=status)

    def send_problem(self, exc, status, *, title=None):
        detail = f"{type(exc).__name__}: {exc}"
        raw = str(exc).casefold()
        if title is None:
            title = str(exc) if status < 500 else "操作未能完成；请查看技术详情。"
        if isinstance(exc, BlockingIOError):
            next_steps = ["等待当前运行结束后刷新状态，再重试此操作。"]
        elif isinstance(exc, FileNotFoundError) or any(token in raw for token in ("path", "路径", "file not found", "does not exist")):
            next_steps = ["检查文件或运行目录是否仍存在，并确认输入位于当前工作目录允许的范围内。"]
        elif isinstance(exc, ModuleNotFoundError) or "no module named" in raw or "missing package" in raw:
            next_steps = ["根据技术详情安装或切换到包含所需依赖的 Python 环境，然后重启工作台。"]
        elif any(token in raw for token in ("count", "计数", "integer", "整数", "negative", "非负", "fractional", "raw")):
            next_steps = ["确认输入为原始、有限、非负整数计数；若数据经过归一化，请改用原始计数矩阵或对应 counts layer。"]
        elif any(token in raw for token in ("condition", "条件", "donor", "供体", "group", "分组", "replicate", "重复")):
            next_steps = ["重新核对条件列、条件值和供体映射，并按实验设计确认每组独立重复数。"]
        else:
            next_steps = ["核对输入与技术详情；若能稳定复现，请保留运行目录和错误信息供进一步排查。"]
        return self.send_json({"error": title, "detail": detail, "next_steps": next_steps}, status)

    def dispatch(self, method):
        # Switching roots and invalidating old page tokens is one operation.
        # No request may resolve a run or input against mixed directory roots.
        with self.server.service.directory_guard:
            self.dispatch_locked(method)

    def dispatch_locked(self, method):
        try:
            if self.headers.get("Host") not in self.server.allowed_hosts:
                return self.send_json({"error": "Invalid Host"}, 403)
            origin = self.headers.get("Origin")
            if origin and origin not in {"http://" + host for host in self.server.allowed_hosts}:
                return self.send_json({"error": "Cross-origin requests are not allowed"}, 403)
            if self.headers.get("Sec-Fetch-Site") == "cross-site":
                return self.send_json({"error": "Cross-site requests are not allowed"}, 403)
            url = urlsplit(self.path)
            path = url.path
            if not path.startswith("/api/"):
                if method != "GET":
                    return self.send_json({"error": "Method not allowed"}, 405)
                assets = {"/": ("index.html", "text/html; charset=utf-8"),
                          "/app.js": ("app.js", "text/javascript; charset=utf-8"),
                          "/results.js": ("results.js", "text/javascript; charset=utf-8"),
                          "/style.css": ("style.css", "text/css; charset=utf-8")}
                if path not in assets:
                    return self.send_json({"error": "Not found"}, 404)
                name, content_type = assets[path]
                body = files("eacbp.webui").joinpath("static", name).read_bytes()
                if name == "index.html":
                    body = body.replace(b"__EACBP_TOKEN__", self.server.token.encode("ascii"))
                return self.send_content(body, content_type)
            if not secrets.compare_digest(self.headers.get("X-EACBP-Token", ""), self.server.token):
                return self.send_json({"error": "页面会话已过期，请刷新浏览器"}, 403)
            query = parse_qs(url.query)
            service = self.server.service
            if method == "GET":
                if path == "/api/settings":
                    return self.send_json(service.settings())
                if path == "/api/files":
                    return self.send_json(service.browse(query.get("path", ["."])[0]))
                if path == "/api/runs":
                    return self.send_json({"runs": service.manager.runs(), "jobs": service.manager.jobs()})
                if path == "/api/run":
                    return self.send_json(service.manager.detail(query.get("id", [""])[0]))
                if path == "/api/reuse":
                    return self.send_json(service.reuse_configuration(query.get("id", [""])[0]))
                if path == "/api/results":
                    return self.send_json(service.results(query.get("id", [""])[0]))
                if path == "/api/results/deg":
                    def one(name, default=None):
                        values = query.get(name)
                        if values is None:
                            return default
                        if len(values) != 1:
                            raise ValueError(f"{name} 只能提供一次")
                        return values[0]
                    page_text, size_text = one("page", "1"), one("size", "50")
                    try:
                        page, size = int(page_text), int(size_text)
                    except (TypeError, ValueError) as exc:
                        raise ValueError("page 和 size 必须是整数") from exc
                    threshold = one("fdr_max")
                    try:
                        threshold = float(threshold) if threshold is not None else None
                    except ValueError as exc:
                        raise ValueError("FDR 阈值必须是数字") from exc
                    fold_change = one("abs_log2fc_min")
                    try:
                        fold_change = float(fold_change) if fold_change is not None else None
                    except ValueError as exc:
                        raise ValueError("绝对 log2 fold change 阈值必须是数字") from exc
                    only = one("significant_only", "false")
                    if only not in {"true", "false"}:
                        raise ValueError("significant_only 必须是 true 或 false")
                    return self.send_json(service.deg_table(
                        one("id", ""), task_id=one("task_id", ""), page=page, size=size,
                        query=one("q", ""), fdr_max=threshold, significant_only=only == "true",
                        abs_log2fc_min=fold_change))
                if path == "/api/results/deg-tables":
                    return self.send_json(service.deg_tables(query.get("id", [""])[0]))
                if path == "/api/results/deg.csv":
                    target, name = service.deg_table_csv(query.get("id", [""])[0],
                                                          task_id=query.get("task_id", [""])[0])
                    return self.send_path(target, "text/csv; charset=utf-8", name)
                if path == "/api/export":
                    target, name = service.export_bundle(query.get("id", [""])[0])
                    return self.send_path(target, "application/zip", name)
                if path == "/api/snapshot":
                    target, name = service.snapshot_download(query.get("id", [""])[0])
                    return self.send_path(target, "application/json; charset=utf-8", name)
                if path == "/api/methods":
                    run_id = query.get("id", [""])[0]
                    body = service.methods_report(run_id).encode("utf-8")
                    return self.send_content(body, "text/markdown; charset=utf-8",
                                             download=f"{Path(run_id).name}-methods.md")
                if path == "/api/download":
                    name = query.get("name", [""])[0]
                    if name not in {"report.md", "manifest.json", "config.json", "summary.json"}:
                        raise ValueError("Unsupported download")
                    run = service.manager.run_path(query.get("id", [""])[0])
                    target = run / name
                    if not target.resolve().is_relative_to(run):
                        raise ValueError("Invalid download path")
                    if target.stat().st_size > 32_000_000:
                        raise ValueError("文件超过浏览器下载限制，请从运行目录读取")
                    return self.send_content(target.read_bytes(), "text/plain; charset=utf-8", download=name)
            elif method == "POST":
                if self.headers.get("Content-Type", "").split(";")[0] != "application/json":
                    raise ValueError("Content-Type must be application/json")
                size = int(self.headers.get("Content-Length", "0"))
                if not 0 < size <= 1_000_000:
                    raise ValueError("请求大小必须在 1 MB 以内")
                self.connection.settimeout(15)
                payload = json.loads(self.rfile.read(size))
                if not isinstance(payload, dict):
                    raise ValueError("Expected JSON object")
                if path == "/api/settings":
                    result = service.update_directories(payload)
                    self.server.token = secrets.token_urlsafe(32)
                    return self.send_json(result)
                if path == "/api/dataset":
                    return self.send_json(service.dataset(payload.get("data", "")))
                if path == "/api/design":
                    return self.send_json(service.design(payload))
                if path == "/api/preview":
                    return self.send_json(service.preview(payload))
                if path == "/api/import/preview":
                    return self.send_json(service.preview_import(payload))
                if path == "/api/import/confirm":
                    return self.send_json(service.confirm_import(payload))
                if path == "/api/run":
                    return self.send_json(service.start(payload), 202)
                if path in {"/api/resume", "/api/report"}:
                    run = service.manager.run_path(payload.get("run_id", ""))
                    job = service.manager.submit(path.rsplit("/", 1)[1], run_dir=run)
                    return self.send_json({"job": job}, 202)
            self.send_json({"error": "Not found"}, 404)
        except BlockingIOError as exc:
            self.send_problem(exc, 409)
        except FileNotFoundError as exc:
            self.send_problem(exc, 404)
        except (ValueError, TypeError, KeyError) as exc:
            self.send_problem(exc, 400)
        except Exception as exc:
            self.send_problem(exc, 500)

    def do_GET(self):
        self.dispatch("GET")

    def do_POST(self):
        self.dispatch("POST")


def serve(*, workspace=None, runs_dir=None, port=8765, open_browser=True, settings_file=None):
    from .preferences import startup_directories

    workspace, runs_dir, settings_file = startup_directories(
        workspace=workspace, runs_dir=runs_dir, settings_file=settings_file)
    server = WebServer(WebService(workspace, runs_dir, settings_file=settings_file), port=port)
    url = f"http://127.0.0.1:{server.server_port}/"
    print(f"EACBP WebUI: {url}\nWorkspace: {workspace}\nRuns: {runs_dir}\nCtrl+C stops the UI; background jobs continue.", flush=True)
    if open_browser:
        webbrowser.open(url)
    try:
        server.serve_forever(poll_interval=0.5)
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0
