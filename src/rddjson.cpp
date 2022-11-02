#include <iostream>
#include "../thirdparty/nlohmann/json.hpp"

using json = nlohmann::json;

/*
int main()
{

    // jdEmployees
    json jdEmployees =
    {
        {"firstName","Sean"},
        {"lastName","Brown"},
        {"StudentID",21453},
        {"Department","Computer Sc."}
    };

    // Access the values
    std::string fName = jdEmployees.value("firstName", "oops");
    std::string lName = jdEmployees.value("lastName", "oops");
    int sID = jdEmployees.value("StudentID", 0);
    std::string dept = jdEmployees.value("Department", "oops");
   
    // Print the values
    std::cout << "First Name: " << fName << std::endl;
    std::cout << "Last Name: " << lName << std::endl;
    std::cout << "Student ID: " << sID << std::endl;
    std::cout << "Department: " << dept << std::endl;
             
    return 0;
}
*/
