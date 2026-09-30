@echo off
setlocal
cd /d "%~dp0"
python -m streamlit run app.py
if errorlevel 1 pause
endlocal
