.PHONY : all html latex
all : html latex
latex : html
	cd latex; make ; cp refman.pdf ../Sky3D_doxygen.pdf; cd ..
html :
	doxygen Doxyfile
clean : 
	-@rm -rf latex html *.tmp *_log.txt *.pdf 2>/dev/null || true 
