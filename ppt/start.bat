@echo off
cd /d "%~dp0"
echo.
echo  PPT 已启动，请在浏览器中访问：
echo  http://localhost:8000/index.html
echo.
echo  按 Ctrl+C 关闭服务器
echo.
python -m http.server 8000
