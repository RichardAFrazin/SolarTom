TOMROOT_HEADER = tomroot.h

#CC = clang 
 CC = gcc    # for Linux systems

DEFINES =
 CFLAGS = -Wall -O3 -I/data1/tomography/wcstools-3.7.0/libwcs -L/data1/tomography/wcstools-3.7.0/libwcs

#################################-I/usr/include/sys

LIBS = -lm  -lwcs

MEX = mex
MEX_FLAGS = -argcheck -g -O -outdir ./matlab/

INDENT = astyle
INDENT_FLAGS = --style=kr --indent=spaces=2


#####################################################################

SRC_BUILD = get_orbit.c r3misc.c rots.c grids.c
OBJ_BUILD = $(SRC_BUILD:%.c=%.o)

SRC_SOLVE = fess_hu.c cg_quad.c sparse.c normcalc.c 
OBJ_SOLVE = $(SRC_SOLVE:%.c=%.o)

SRC_AUX = amoeba.c fminbr.c builda.c datetest.c callsolve.c row_extract.c \
          calculate_y.c row_to_col.c fitstest.c \
	  fitstest.c solve.c build_subA.c get_build_opts.c

BIN = builda callsolve_cg callsolve_fess auto_cv_brent compare

all: tomroot.h $(BIN)

tomroot.h:  $(TOMROOT_HEADER) 
	cp $(TOMROOT_HEADER) tomroot.h

builda: $(OBJ_BUILD) builda.o build_subA.o rcs_llist.o llist.o

compare: $(OBJ_BUILD) compare.o

datetest: $(OBJ_BUILD) datetest.o

## print_grid: $(CC) buildA_params.h grids.c print_grid.c -o print_grid

print_grid: $(OBJ_BUILD) print_grid.o

solve_cg.o: headers.h
solve_cg.o: CFLAGS += -UFESSMIN -DCONJGRAD
solve_cg.o: solve.c
	$(CC) -c $(CFLAGS) $(DEFINES) $< -o $@

callsolve_cg: $(OBJ_SOLVE) callsolve.o solve_cg.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

solve_fess.o: headers.h
solve_fess.o: CFLAGS += -UCONJGRAD -DFESSMIN
solve_fess.o: solve.c
	$(CC) -c $(CFLAGS) $(DEFINES) $< -o $@

callsolve_fess: $(OBJ_SOLVE) callsolve.o solve_fess.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

###auto_cv: $(OBJ_SOLVE) auto_cv.o amoeba.o solve_cg.o
###	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

auto_cv_brent: $(OBJ_SOLVE) auto_cv_brent.o fminbr.o cvcalc.o solve_cg.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

cv_brent_fixed: $(OBJ_SOLVE) cv_brent_fixed.o fminbr.o cvcalc.o solve_cg.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

row_extract: row_extract.o sparse.o

calculate_y: calculate_y.o sparse.o


extract_block: extract_block.o sparse.o

llist_test: llist_test.o llist.o

rcs_llist_test: rcs_llist_test.o rcs_llist.o llist.o

#./matlab/get_build_opts.mexglx: get_build_opts.c
#	$(MEX) $(MEX_FLAGS) $<

#./matlab/get_auto_median_opts.mexglx: get_auto_median_opts.c
#	$(MEX) $(MEX_FLAGS) $<

indent:
	$(INDENT) $(INDENT_FLAGS) *.c *.h

%.o : %.c
	$(CC) -c $(CFLAGS) $(DEFINES) $<

% : %.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

clean:
	rm -f *.o $(BIN)

# DO NOT DELETE
headers.h: buildA_params.h solve_cv_params.h  tomroot.h
get_orbit.o: headers.h  buildA_params.h solve_cv_params.h
r3misc.o: headers.h  buildA_params.h solve_cv_params.h
rots.o: headers.h  buildA_params.h solve_cv_params.h
fess_hu.o: headers.h  buildA_params.h solve_cv_params.h
cg_quad.o: headers.h  buildA_params.h solve_cv_params.h
sparse.o: headers.h  buildA_params.h solve_cv_params.h
cvcalc.o: headers.h  buildA_params.h solve_cv_params.h
normcalc.o: headers.h  buildA_params.h solve_cv_params.h
llist.o: llist.h
rcs_llist.o: rcs_llist.h llist.h
builda.o: headers.h  buildA_params.h solve_cv_params.h
compare.o: buildrow.c buildA_params.h
build_subA.o: buildrow.c  buildA_params.h solve_cv_params.h
callsolve.o: headers.h  buildA_params.h solve_cv_params.h
row_extract.o: headers.h  buildA_params.h solve_cv_params.h
solve.o: headers.h buildA_params.h solve_cv_params.h
fminbr.o: headers.h buildA_params.h solve_cv_params.h
cv_brent_fixed.o: headers.h buildA_params.h solve_cv_params.h
auto_cv_brent.o: headers.h buildA_params.h solve_cv_params.h
extract_block.o: headers.h buildA_params.h solve_cv_params.h 
print_grid.o: buildA_params.h grids.c print_grid.c

