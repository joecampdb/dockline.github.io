@echo off
echo Starting DiffDock for imatinib_Abl...
cd /d "C:\Users\ltjjp\OneDrive\Desktop\joecam\diffdock-test\DiffDock"
echo Current directory: %CD%
"C:\Users\ltjjp\anaconda3\envs\diffdock\python.exe" -m inference --protein_path "C:\Users\ltjjp\OneDrive\Desktop\joecam\diffdock-test\pdb_clean\imatinib_Abl_clean.pdb" --ligand "C:\Users\ltjjp\OneDrive\Desktop\joecam\diffdock-test\imatinib.sdf" --out_dir "C:\Users\ltjjp\OneDrive\Desktop\joecam\diffdock-test\results\imatinib_Abl" --samples_per_complex 10 --actual_steps 18 --no_final_step_noise
echo Done with exit code: %ERRORLEVEL%
