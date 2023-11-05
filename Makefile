NC = nvcc

SRC1g2m = cuda-1g-wsm2.cu

OUT1g2m = cuda-1g-wsm2.out


2m:$(SRC1g2m)
	$(NC) $(SRC1g2m) -o $(OUT1g2m)

2mrun:2m
	./$(OUT1g2m)

clean:
	rm $(OUT1g2m)
