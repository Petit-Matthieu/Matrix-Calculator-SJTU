@echo off
cd /d "%~dp0"
echo.
echo  Server started, open in browser:
echo  http://localhost:8000/index.html
echo.
echo  Press Ctrl+C to stop
echo.
python -m http.server 8000
