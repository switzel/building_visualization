all: char0_vertex.obj char0_chamber.obj char2_vertex.obj char2_chamber.obj

clean:
	rm char0_vertex.obj
	rm char0_chamber.obj
	rm char2_vertex.obj
	rm char2_chamber.obj

%.obj: %.simp
	python3 augment_embedding.py --transform --wavefront $< $@
