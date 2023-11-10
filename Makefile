NC = nvcc

CFLAGS = 

SRC1g2m = cuda-1g-wsm2.cu

OBJ_1g2m = cuda-1g-wsm2.o

SRC_SSSP = sssp.cu

OUT_SSSP = sssp.out


sssp:$(OBJ_1g2m) $(SRC_SSSP) 
	$(NC) $(CFLAGS) $^ -o $(OUT_SSSP)

$(OBJ_1g2m):$(SRC1g2m)
	$(NC) $(CFLAGS) -c $^ -o $@

run:sssp
	./$(OUT_SSSP)

clean:
	rm $(OUT_SSSP) $(OBJ_1g2m)
