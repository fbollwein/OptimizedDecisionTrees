CC=g++
CC_FLAGS = -std=c++17 -O3 -c -fpic -Wno-ignored-attributes
LTO=-flto
SWIGV= /usr/bin/swig

MKDIR_P = mkdir -p
BUILD_DIR = build
LIB_DIR = lib

BOOST_DIR=
EIGEN_DIR=
LIBLINEAR_DIR=
CLP_DIR=
CTPL_DIR=
EXPRTK_DIR=
EIGENMVN_DIR=

USEGUROBI =
GUROBI_INC=
GUROBI_LIB= 

ifneq ($(EIGEN_DIR),)
	EIGEN_INC=-I$(EIGEN_DIR)
else
	EIGEN_INC=
endif

ifneq ($(BOOST_DIR),)
	BOOST_INC=-I$(BOOST_DIR)
else
	BOOST_INC=
endif

ifneq ($(CLP_DIR),)
	CLP_INC=-I$(CLP_DIR)/include/coin
	CLPLIB=-L$(CLP_DIR)/lib -lClpSolver -lClp -lCoinUtils -Wl,-rpath,$(CLP_DIR)/lib
else
	CLP_INC=-I$(CLP_DIR)/include/coin
	CLPLIB=-L$(CLP_DIR)/lib -lClpSolver -lClp -lCoinUtils -Wl,-rpath,$(CLP_DIR)/lib
endif

ifneq ($(CTPL_DIR),)
	CTPL_INC=-I$(CTPL_DIR)
else
	CTPL_INC=
endif

ifneq ($(EXPRTK_DIR),)
	EXPRTK_INC=-I$(EXPRTK_DIR)
else
	EXPRTK_INC=
endif

ifneq ($(LIBLINEAR_DIR),)
	LIBLINEAR_INC=-I$(LIBLINEAR_DIR)
else
	LIBLINEAR_INC=
endif

ifneq ($(EIGENMVN_DIR),)
	EIGENMVN_INC=-I$(EIGENMVN_DIR)
else
	EIGENMVN_INC=
endif

SWIGINCLUDE= -I./
PY_CFLAGS  = $(shell python3-config --includes)

INCLUDE= $(EIGEN_INC) $(CLP_INC) $(CTPL_INC) $(EXPRTK_INC) $(LIBLINEAR_INC) $(EIGENMVN_INC) $(GUROBI_INC) $(BOOST_INC)

TREE=DecisionTree
SO = so
SWIG = $(BUILD_DIR)/$(TREE)_wrap.o
FOREST_SWIG = $(BUILD_DIR)/$(FOREST)_wrap.o
SOURCES = $(wildcard src/*.cpp)
HEADERS = $(wildcard src/*.h)
OBJECTSTMP = $(patsubst src/%,$(BUILD_DIR)/%,$(SOURCES))
OBJECTS = $(OBJECTSTMP:.cpp=.o)

LIBLINEAR = $(LIBLINEAR_DIR)/linear.o
LIBNEWTON= $(LIBLINEAR_DIR)/newton.o
LINEARLIBS= $(LIBLINEAR_DIR)/blas/blas.a

.PHONY: directories all clean
 
$(SO): $(OBJECTS) $(SWIG) $(LIBLINEAR) $(LIBNEWTON) $(LINEARLIBS)
	$(CC) $(LTO) -shared $(OBJECTS) $(LIBLINEAR) $(LIBNEWTON) $(SWIG) -o lib/_$(TREE).so $(GUROBI_LIB) $(CLPLIB) -lm $(LINEARLIBS)

$(SWIG):
	$(SWIGV) -c++ -python -outdir lib -o src/$(TREE)_wrap.cxx swig/$(TREE).i
	$(CC) $(CC_FLAGS) $(LTO) src/$(TREE)_wrap.cxx -o $(SWIG) $(PY_CFLAGS) $(INCLUDE) $(SWIGINCLUDE)

$(BUILD_DIR)/%.o: src/%.cpp $(LIBLINEAR) $(LIBNEWTON) $(LINEARLIBS) | directories
	$(CC) $(CC_FLAGS) $(LTO) $(INCLUDE) $< -o $@ $(USEGUROBI)

$(LIBLINEAR) $(LIBNEWTON) $(LINEARLIBS):
	make -C $(LIBLINEAR_DIR)
 
directories: $(BUILD_DIR) $(LIB_DIR)

$(BUILD_DIR):
	${MKDIR_P} $(BUILD_DIR)

$(LIB_DIR):
	${MKDIR_P} $(LIB_DIR)

clean:
	rm -f build/*.o && rm -f *.h.gch && rm -f src/*wrap.cxx && rm -f lib/*.so
	rm -rf $(BUILD_DIR)
	rm -rf $(LIB_DIR)
