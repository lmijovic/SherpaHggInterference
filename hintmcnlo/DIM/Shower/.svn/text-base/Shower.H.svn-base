#ifndef DIM__Shower__Shower_H
#define DIM__Shower__Shower_H

#include "DIM/Tools/Parton.H"
#include "DIM/Shower/Kernel.H"
#include "DIM/Tools/Amplitude.H"

namespace ATOOLS {
  class Data_Reader;
}
namespace MODEL {
  class Model_Base;
}
namespace PDF {
  class PDF_Base;
  class ISR_Handler;
}

namespace DIM {

  class Gamma;

  class Shower {
  public:

    typedef std::map<ATOOLS::Flavour,Kernel_Vector> SKernel_Map;

    typedef std::map<ATOOLS::Flavour,Kernel*> EKernel_Map;
    typedef std::map<ATOOLS::Flavour,EKernel_Map> SEKernel_Map;
    typedef std::map<int,SEKernel_Map> Kernel_Map;

  private:

    MODEL::Model_Base *p_model;
    PDF::PDF_Base     *p_pdf[2];

    Gamma *p_gamma;

    Kernel_Vector m_cks;
    SKernel_Map   m_sks;
    Kernel_Map    m_kmap;

    double m_tmin[2], m_cplfac[2], m_rsf;
    double m_weight, m_oef, m_pdfmin;

    int m_kfac, m_cpl;

    unsigned int m_maxem;

    Splitting m_s;

    void AddKernel(Kernel *const k);

    Splitting GeneratePoint(Parton &p,const double &t);
    Splitting GeneratePoint(const Amplitude &a,const double &t);

  public:

    Shower();

    ~Shower();

    bool Init(MODEL::Model_Base *const model,
	      PDF::ISR_Handler *const isr,
	      ATOOLS::Data_Reader *const read);

    void SetMS(ATOOLS::Mass_Selector *const ms);

    int Evolve(Amplitude &a,unsigned int &nem);

    double GetXPDF(const double &x,const double &Q2,
		   const ATOOLS::Flavour &fl,const int b) const;

    Kernel *GetKernel(const Splitting &s,const int mode) const;

    inline MODEL::Model_Base *Model() const { return p_model; }

    inline double TMin(const int i) const { return m_tmin[i]; }

    inline double CplFac(const int i) const { return m_cplfac[i]; }

    inline double PDFMin() const { return m_pdfmin; }

    inline double GetWeight() const { return m_weight; }

    inline int KFactorScheme() const  { return m_kfac; }
    inline int CouplingScheme() const { return m_cpl;  }

    inline double MuRFactor() const { return m_rsf; }

    inline const Splitting &LastSplitting() const { return m_s; }

    inline void SetGamma(Gamma *const gamma) { p_gamma=gamma; }

    inline double OEF() const { return m_oef; }

  };//end of class Shower

}// end of namespace DIM

#endif
