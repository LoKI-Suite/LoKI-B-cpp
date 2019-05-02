#include <iostream>
#include <string>
#include <lua/lua.hpp>
#include <LuaBridge/LuaBridge.h>
#include <LuaBridge/Vector.h>

using namespace luabridge;

double electronCharge = 1.6021766208e-19;
double boltzmann = 1.38064852e-23;

void print_error(lua_State* state)
{
    // The error message is at the top of the stack.
    // Fetch it, print it and then pop it off the stack.
    std::cout << lua_tostring(state, -1) << std::endl;
    lua_pop(state, 1);
}

void construct_argument_vector(std::vector<double> &arguments, int argc, char ** argv)
{
    // loop through all arguments (after the first, which is the function name)
    for (uint8_t i = 2; i < argc; ++i)
    {
        // cast the argument to a double and add it to the argument vector
        arguments.emplace_back(std::strtod(argv[i], nullptr));
    }
}

int main(int argc, char ** argv)
{
    // check if an argument is provided
    if (argc < 2)
    {
        std::cout << "Incorrect number of arguments. At least one argument expected." << std::endl;
        return 1;
    }

    // build path string
    std::string lua_file = "../lua_test/scripts/";
    lua_file += argv[1];
    lua_file += ".lua";

    // create a new lua state
    lua_State * state = luaL_newstate();

    // make standard libs available to newly constructed state
    luaL_openlibs(state);

    // add a namespace "PhysConst" containing constants to be used in lua scripts
    getGlobalNamespace(state)
            .beginNamespace("PhysConst")
             .addVariable("boltzmann", &boltzmann, false)
             .addVariable("electronCharge", &electronCharge, false)
            .endNamespace();

    // load the file and clean the stack (?)
    // if the file could not be opened then print the error and return
    if (luaL_loadfile(state, lua_file.c_str()) || lua_pcall(state, 0, 0, 0))
    {
        print_error(state);
        return 2;
    }

    // make a vector containing the function arguments
    std::vector<double> arguments;

    // convert the char* array (argv) to a vector of doubles
    construct_argument_vector(arguments, argc, argv);

    // obtain the function from the file (its title is stored in the first argument)
    LuaRef func = getGlobal(state, argv[1]);

    // call the function and cast its result to double
    auto result = func(arguments).cast<double>();

    // print the result
    std::cout << "result: " << result << std::endl;

    return 0;
}