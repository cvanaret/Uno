for f in *.mod; do
	filename=$(basename "$f")
	filename="${filename%.*}"
	/media/data/Bureau/NEOS_solvers/ampl -og$filename $f;
done
