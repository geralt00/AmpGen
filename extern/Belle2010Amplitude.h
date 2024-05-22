#ifndef __Belle2010Amplitude__
#define __Belle2010Amplitude__

#include <string>
#include "extern/DtoKpipiAmplitude.h"

class Belle2010Amplitude : public DtoKpipiAmplitude {

    public:
        Belle2010Amplitude( std::string kaon_type, std::string cp_mult_dir);
        ~Belle2010Amplitude(); // non-virtual: should not be a base class

        virtual dcomplex get_amplitude(double sKp, double sKm) const;
        virtual std::string get_amp_prefix() const {return "Belle2010_";};

    protected:
        void set_amplitude_factors();
        double _r_value[8];
        double _delta[8];


};


#endif