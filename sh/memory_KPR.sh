# Set your other parameters here
PROJECT="memory_response"
SUBPROJECT="Z_NA_dynamics"
PYTHON_SCRIPT="/Users/robertomorantovar/Documents/GitHub/Immune_System/Codes/python/_general_simulations/MeRiM.py"

for p in $(seq 1.0 0.5 4.0); do
    # Run primary response
    python "$PYTHON_SCRIPT" --p $p --pmem 1.0 --secondary 0 --pro $PROJECT --subpro $SUBPROJECT

    # Find output directory and file (adjust if needed)
    OUTDIR=$(find /Users/robertomorantovar/Dropbox/Research/Immune_system/$PROJECT/$SUBPROJECT -type d -name "L0-*" | sort | tail -n 1)
    ACTIVATED_FILE=$(find "$OUTDIR" -type f -name "activated_repertoire.csv" | head -n 1)

    # Run secondary responses for each pmem
    for pmem in $(seq 1.0 0.5 4.0); do
        python "$PYTHON_SCRIPT" --p $p --pmem $pmem --secondary 1 --input_memory_file "$ACTIVATED_FILE" --pro $PROJECT --subpro $SUBPROJECT
    done
done