# 定义编译器和编译选项
CC = g++
CFLAGS = -fPIC -g -O3 -std=c++17 -fopenmp

# 定义源文件和对象文件
SRCS_NUPBO = ls_assigner/NuPBO/basis_pms.cpp \
             ls_assigner/NuPBO/heuristic.cpp \
             ls_assigner/NuPBO/init.cpp \
             ls_assigner/NuPBO/MaxSATFormula.cpp \
             ls_assigner/NuPBO/my_class.cpp \
             ls_assigner/NuPBO/nuPBOFunction.cpp \
             ls_assigner/NuPBO/parse_arguments.cpp \
             ls_assigner/NuPBO/settings.cpp 

SRCS_LSAssigner = ls_assigner/buildmodel/model.cpp \
                  ls_assigner/LSAssigner.cpp 

OBJS_NUPBO = $(SRCS_NUPBO:.cpp=.o)
OBJS_LSAssigner = $(SRCS_LSAssigner:.cpp=.o)

# 定义动态库文件
LIBLSAssigner = liblsassigner.so

# 定义静态库文件
LIBNUPBO = libnupbo.a

# 定义编译规则
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# 定义生成动态库规则
$(LIBLSAssigner): $(OBJS_LSAssigner) $(LIBNUPBO)
	$(CC) $(CFLAGS) -shared -o $(LIBLSAssigner) $(OBJS_LSAssigner) -L. -lnupbo

# 定义生成静态库规则
$(LIBNUPBO): $(OBJS_NUPBO)
	ar rcs $(LIBNUPBO) $(OBJS_NUPBO)

# 默认目标，编译动态库和静态库
all: $(LIBLSAssigner) $(LIBNUPBO)

# 定义清理规则
clean:
	rm -f $(OBJS_LSAssigner) $(OBJS_NUPBO) $(LIBLSAssigner) $(LIBNUPBO)

# make
# g++ -fPIC -g -O3 -std=c++17 -fopenmp -c test.cpp -o test.o
# g++ -fopenmp -m64 -g -O3 -std=c++17 -o TA main.cpp test.o -L. -llsassigner
# ./TA ispd18_test5_metal5.input.def.ta.txt  >out.txt