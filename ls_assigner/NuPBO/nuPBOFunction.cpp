#include "basis_pms.h"

#include "cmdline.h"
#include "parse_arguments.h"
#include "settings.h"

#include "my_class.h"
#include <iostream>
#include <sstream>
#include <string>

// extern long long total_step;
namespace lsa
{
void interrupt(int sig,
               externvariables &obj) //中断信号触发，输出最佳解，终止程序
{
    print_best_solution(obj);
    // free_memory();
    exit(10);
}

int nuPBOFunction(int argc, char *argv[], std::string outputFileName)
{
    //处理全局变量
    externvariables externVars;
    externVars.outputFile.open(outputFileName);
    if (!externVars.outputFile.is_open())
        {
            // 处理文件打开失败的情况
            std::cerr << "Failed to open file: " << outputFileName << std::endl;
        }

    default_algorithm_settings(externVars); //默认设置
    Config cfg;                             // Config结构体，其中包含输入文件的路径，种子值（默认为1）
    if (!parse_args(argc, argv, cfg))       //把两个参数储存在结构体cfg中
        {
            std::cout << "parse args failed." << std::endl;
            return -1;
        }

    // srand((unsigned)time(NULL));
    srand(cfg.seed);          //生成随机序列
    start_timing(externVars); //起始时间存在&start_time里

    // signal(SIGTERM, interrupt);//检测到SIGTERM信号则输出最佳解并终止

    externVars.is_print = 0; //求解的时候用到
    std::vector<int> init_solution;
    build_instance(cfg.input_file.c_str(), externVars); //通过文件名构建实例
    settings(externVars);                               //一些变量的默认值
    parse_arguments(argc, argv, externVars);            //可以通过参数修改变量的值
    local_search(init_solution, cfg.input_file.c_str(),
                 externVars); //给定初始解，通过局部搜索获得更好的解

    // s.simple_print();
    print_best_solution(externVars);
    free_memory(externVars);

    externVars.outputFile.close();

    return 0;
}
} // namespace lsa