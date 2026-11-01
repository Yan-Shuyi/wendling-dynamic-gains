@echo off
REM ============================================================
REM build_figures.bat - recompile the two TikZ figures
REM
REM   Fig. 1   fig_model_structure.tex   Wendling microcircuit schematic
REM   Fig. 2   fig_classifier_tree.tex   EEG state decision tree
REM
REM   The MATLAB figures (Fig. 3 to Fig. 8) are NOT built here;
REM   run run_all, or the individual run_*.m scripts, in code\.
REM ============================================================
cd /d "%~dp0"

echo [1/2] fig_model_structure.tex
pdflatex -interaction=nonstopmode -halt-on-error fig_model_structure.tex >nul
if errorlevel 1 goto fail
echo       OK

echo [2/2] fig_classifier_tree.tex
pdflatex -interaction=nonstopmode -halt-on-error fig_classifier_tree.tex >nul
if errorlevel 1 goto fail
echo       OK

del /q *.aux >nul 2>&1
del /q *.log >nul 2>&1

echo.
echo DONE  -^>  fig_model_structure.pdf , fig_classifier_tree.pdf
goto end

:fail
echo.
echo FAILED - open the .log file in this folder to see the error

:end
