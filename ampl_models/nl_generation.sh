for f in *.mod; do
	filename=$(basename "$f")
	filename="${filename%.*}"
	~/Desktop/ampl.linux64/ampl -og$filename $f;
done
