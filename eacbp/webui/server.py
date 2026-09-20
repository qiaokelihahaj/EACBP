"""Loopback-only HTTP server; no additional web framework or frontend build."""
from __future__ import annotations

from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib.resources import files
import json
from pathlib import Path
import secrets
from urllib.parse import parse_qs, urlsplit
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

    def send_json(self, value, status=200):
        self.send_content(json.dumps(value, ensure_ascii=False, allow_nan=False).encode("utf-8"), status=status)

    def dispatch(self, method):
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
                if path == "/api/dataset":
                    return self.send_json(service.dataset(payload.get("data", "")))
                if path == "/api/preview":
                    return self.send_json(service.preview(payload))
                if path == "/api/run":
                    return self.send_json(service.start(payload), 202)
                if path in {"/api/resume", "/api/report"}:
                    run = service.manager.run_path(payload.get("run_id", ""))
                    job = service.manager.submit(path.rsplit("/", 1)[1], run_dir=run)
                    return self.send_json({"job": job}, 202)
            self.send_json({"error": "Not found"}, 404)
        except BlockingIOError as exc:
            self.send_json({"error": str(exc)}, 409)
        except FileNotFoundError as exc:
            self.send_json({"error": str(exc)}, 404)
        except (ValueError, TypeError, KeyError) as exc:
            self.send_json({"error": str(exc)}, 400)
        except Exception as exc:
            self.send_json({"error": f"{type(exc).__name__}: {exc}"}, 500)

    def do_GET(self):
        self.dispatch("GET")

    def do_POST(self):
        self.dispatch("POST")


def serve(*, workspace=None, runs_dir=None, port=8765, open_browser=True):
    workspace = Path(workspace or Path.cwd()).expanduser().resolve()
    if not workspace.is_dir():
        raise ValueError("工作目录不存在")
    runs_dir = Path(runs_dir).expanduser().resolve() if runs_dir else workspace / "outputs" / "runs"
    server = WebServer(WebService(workspace, runs_dir), port=port)
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
