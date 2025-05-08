#include <iostream>
#include <string>
#include <vector>

// define a struct that represents a row of data
struct DataRow {
    std::string str;
    float f1;
    int i1;
};

int main() {
    // create a vector to store the data rows
    std::vector<DataRow> data;

    // loop over a function and produce new rows of data
    for (int i = 0; i < N; i++) {
        // create a new data row and add it to the vector
        DataRow row;
        row.str = "abc";
        row.f1 = 1.23;
        row.i1 = 456;
        data.push_back(row);
    }

    // print the data rows to a binary file
    std::ofstream out("data.bin", std::ios::binary);
    for (const auto& row : data) {
        // write the string to the file
        out.write(row.str.c_str(), row.str.size());

        // write the float and int to the file
        out.write(reinterpret_cast<const char*>(&row.f1), sizeof(row.f1));
        out.write(reinterpret_cast<const char*>(&row.i1), sizeof(row.i1));
    }

    return 0;
}
