@echo off
cd /d "C:\Users\ltjjp\OneDrive\Desktop\joecam\diffdock-test\DiffDock"
"C:\Users\ltjjp\anaconda3\envs\diffdock\python.exe" -m inference --config default_inference_args.yaml --protein_ligand_csv "C:\Users\ltjjp\OneDrive\Desktop\joecam\diffdock-test\input_remaining_docking.csv" --out_dir "C:\Users\ltjjp\OneDrive\Desktop\joecam\diffdock-test\results" --samples_per_complex 10 --batch_size 10 --actual_steps 18 --no_final_step_noise
