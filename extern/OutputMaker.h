#ifndef __OUTPUT_MAKER__
#define __OUTPUT_MAKER__

#include <string>
#include "DtoKpipiAmplitude.h"

class OutputMaker {
public:
    OutputMaker();
    ~OutputMaker();

    void write_output(const DtoKpipiAmplitude& amp, const std::string& output_dir, const int numpoints) const;

private:
    std::string get_file_name(const DtoKpipiAmplitude& amp, const int numpoints) const;
};

#endif