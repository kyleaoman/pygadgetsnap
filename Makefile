.PHONY : py_wrap
py_wrap: gadget_binary_snap/snapclass.cpp gadget_binary_snap/snapclass.h gadget_binary_snap/snapclass.i Makefile
	swig -c++ -python -o gadget_binary_snap/snap_wrap.cpp gadget_binary_snap/snapclass.i
