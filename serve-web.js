#!/usr/bin/env node

const http = require('node:http');
const fs = require('node:fs');
const path = require('node:path');
const { URL } = require('node:url');

const HOST = process.env.HOST || '127.0.0.1';
const PORT = Number(process.env.PORT || 8000);
const WEB_ROOT = path.resolve(__dirname, 'web');

const CONTENT_TYPES = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'application/javascript; charset=utf-8',
  '.mjs': 'application/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.png': 'image/png',
  '.jpg': 'image/jpeg',
  '.jpeg': 'image/jpeg',
  '.gif': 'image/gif',
  '.svg': 'image/svg+xml',
  '.ico': 'image/x-icon',
  '.txt': 'text/plain; charset=utf-8',
  '.map': 'application/json; charset=utf-8'
};

function contentTypeFor(filePath) {
  return CONTENT_TYPES[path.extname(filePath).toLowerCase()] || 'application/octet-stream';
}

function resolveFromWebRoot(requestUrl) {
  const parsed = new URL(requestUrl, 'http://localhost');
  let pathname = decodeURIComponent(parsed.pathname);

  if (pathname === '/' || pathname === '') {
    pathname = '/index.html';
  }

  const candidate = path.resolve(WEB_ROOT, `.${pathname}`);
  if (!candidate.startsWith(WEB_ROOT)) {
    return null;
  }

  return candidate;
}

function sendText(res, status, text) {
  res.writeHead(status, { 'Content-Type': 'text/plain; charset=utf-8' });
  res.end(text);
}

const server = http.createServer((req, res) => {
  if (!req.url) {
    return sendText(res, 400, 'Bad Request');
  }

  if (req.method !== 'GET' && req.method !== 'HEAD') {
    res.setHeader('Allow', 'GET, HEAD');
    return sendText(res, 405, 'Method Not Allowed');
  }

  const resolved = resolveFromWebRoot(req.url);
  if (!resolved) {
    return sendText(res, 403, 'Forbidden');
  }

  fs.stat(resolved, (err, stats) => {
    if (err) {
      return sendText(res, 404, 'Not Found');
    }

    let filePath = resolved;
    if (stats.isDirectory()) {
      filePath = path.join(resolved, 'index.html');
    }

    fs.stat(filePath, (statErr, fileStats) => {
      if (statErr || !fileStats.isFile()) {
        return sendText(res, 404, 'Not Found');
      }

      const headers = {
        'Content-Type': contentTypeFor(filePath),
        'Content-Length': fileStats.size,
        'Cache-Control': 'no-cache'
      };

      res.writeHead(200, headers);

      if (req.method === 'HEAD') {
        return res.end();
      }

      const stream = fs.createReadStream(filePath);
      stream.on('error', () => sendText(res, 500, 'Internal Server Error'));
      stream.pipe(res);
    });
  });
});

server.listen(PORT, HOST, () => {
  console.log(`Serving ${WEB_ROOT}`);
  console.log(`Open: http://${HOST}:${PORT}`);
});

