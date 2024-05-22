#ifndef __DtoKpipiAmplitude__
#define __DtoKpipiAmplitude__

#include <string>
#include <vector>
#include <complex>


typedef std::complex<double> dcomplex;

class DtoKpipiAmplitude{

    public:

        DtoKpipiAmplitude(std::string kaon_type);
        ~DtoKpipiAmplitude(); // non-virtual: should not be a base class

        virtual dcomplex get_amplitude(double sKp, double sKm) const = 0;

        bool is_valid_point(double x, double y) const;

        double get_min_s() const {return (_mKs + _mPi)*(_mKs + _mPi);};
        double get_max_s() const {return (_mD  - _mPi)*(_mD  - _mPi);};

        std::string get_kaon_type() const {return _kaon_type;};

        virtual std::string get_amp_prefix() const {return "";};

    protected:
        virtual void set_amplitude_factors() = 0;
        dcomplex get_amplitude_factor(double r, double delta) const;

        const std::string _kaon_type;

        // Factors between KS and KL amplitudes, A(KL)/A(KS) for the resonances
        dcomplex _CF_factor; // Cabibbo favoured
        dcomplex _DCS_factor; // Doube Cabibbo suppressed
        dcomplex _nonRes_factor; // Non-resonant part
        std::vector<dcomplex> _pipi_factors; // vector holding the individual pipi-resonance values

        // Particle masses
        const double _mD;
        const double _mKs;
        const double _mPi;

};


#endif