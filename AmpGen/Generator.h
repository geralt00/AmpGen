#ifndef AMPGEN_GENERATOR_H
#define AMPGEN_GENERATOR_H

#include "AmpGen/EventList.h"
#if ENABLE_AVX 
  #include "AmpGen/EventListSIMD.h"
#endif
#include "AmpGen/simd/utils.h"
#include "AmpGen/EventType.h"
#include "AmpGen/PhaseSpace.h"
#include "AmpGen/TreePhaseSpace.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/ProfileClock.h"
#include "AmpGen/ProgressBar.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/MetaUtils.h"

namespace AmpGen
{
  template <typename phaseSpace_t = PhaseSpace>
    class Generator
    {
        #if ENABLE_AVX
          using eventlist_t = EventListSIMD;
        #else
          using eventlist_t = EventList;
        #endif

      private:
        EventType    m_eventType;
        phaseSpace_t m_gps;
        size_t       m_generatorBlock = {5000000};
        TRandom*     m_rnd            = {gRandom};
        bool         m_normalise      = {true};

      public:
        template <typename... ARGS> explicit Generator( const ARGS&... args ) : m_gps(args...)
        {
          m_eventType = m_gps.eventType();
          if( m_rnd != gRandom ) setRandom( m_rnd );
        }
  
        phaseSpace_t phsp() { return m_gps; }
  
        void setRandom( TRandom* rand )
        {
          m_rnd = rand;
          m_gps.setRandom( m_rnd );
        }
        void setBlockSize( const size_t& blockSize ) { m_generatorBlock = blockSize; }
        void setNormFlag( const bool& normSetting )  { m_normalise = normSetting; }

        void fillEventListPhaseSpace( eventlist_t& events, const size_t& N)
        {
          if constexpr( std::is_same<phaseSpace_t, PhaseSpace>::value )
          {
            constexpr auto w = utils::size<real_v>::value; 
            for( unsigned i = 0 ; i != N ; ++i ){
              double* addr = reinterpret_cast<double*>( events.block( i/w  ))+ i % w;
              m_gps.fill(addr, w);     
            }     
          }
          else {
            auto it = events.begin();
            while( it != events.end() )
            {
              *it = m_gps.makeEvent();
              ++it;
            }
          }
        }
        template <typename pdf_t> 
        double getMax(const eventlist_t& events, pdf_t& pdf ) const 
        {
          double max = 0.;
          for ( const auto& evt : events ) 
          {
            auto value             = evt.genPdf();
            if( std::isnan(value) ){ 
              ERROR( "PDF for event is nan: " << value  );
              evt.print(); 
              //pdf.debug( evt ); 
            }
            else if ( value > max ) max = value;
          }
          DEBUG( "Returning normalisation constant = " << max ); 
          return max;
        }

        template <typename eventList_t, typename pdf_t> void fillEventList( pdf_t& pdf, eventList_t& list, const size_t& N)
        {
          if ( m_rnd == nullptr ) 
          {
            ERROR( "Random generator not set!" );
            return;
          }
          double maxProb   = m_normalise ? 0 : 1;
          auto size0       = list.size();
          double totalGenerated = 0; 
          pdf.reset( true );
          ProgressBar pb(60, detail::trimmedString(__PRETTY_FUNCTION__) );
          ProfileClock t_phsp, t_eval, t_acceptReject, t_total;
          std::vector<bool> efficiencyReport(m_generatorBlock, false); 

          while ( list.size() - size0 < N ) {
            eventlist_t mc( m_eventType );
            mc.resize(m_generatorBlock);
            t_phsp.start();
            fillEventListPhaseSpace(mc, m_generatorBlock);
            t_phsp.stop();
            t_eval.start();
            pdf.setEvents(mc);
            pdf.prepare();
            auto previousSize = list.size();
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for ( size_t block=0; block < mc.nBlocks(); ++block )
            {
              mc.setGenPDF(block, pdf(mc.block(block), block) / mc.genPDF(block) );
            }
            maxProb = maxProb == 0 ? 1.2 * getMax(mc, pdf) : maxProb; 
            DEBUG( "Norm: " << maxProb );            
           // if constexpr ( std::is_same<phaseSpace_t, TreePhaseSpace>::value ) m_gps.recalculate_weights(mc); 

            t_eval.stop();
            t_acceptReject.start(); 
            totalGenerated += mc.size();
            for(const auto& event : mc)
            { 
              if ( event.genPdf()  > maxProb ) {
                std::cout << std::endl; 
                WARNING( "PDF value exceeds norm value: " << event.genPdf() << " > " << maxProb );
                event.print();
              }
              if ( event.genPdf() > maxProb * m_rnd->Rndm() ){
                list.push_back(event);
                list.rbegin()->setGenPdf( pdf(event) );
                efficiencyReport[event.index()] = true; 
              }
              else efficiencyReport[event.index()] = false; 
              if ( list.size() - size0 == N ) break; 
            }
            /// if constexpr ( std::is_same<phaseSpace_t, TreePhaseSpace>::value ) maxProb = 0; 
            t_acceptReject.stop(); 
            double efficiency = 100. * ( list.size() - previousSize ) / (double)m_generatorBlock;
            pb.print( double(list.size()) / double(N), " ε[gen] = " + mysprintf("%.4f",efficiency) + "% , " + std::to_string(int(t_total.count()/1000.))  + " seconds" );
            if ( list.size() == previousSize ) {
              ERROR( "No events generated, PDF: " << type_string<pdf_t>() << " is likely to be malformed" );
              break;
            }
          } 
          pb.finish();
          t_total.stop();
          INFO("Generated " << N << " events in " << t_total << " ms");
          INFO("Generating phase space : " << t_phsp         << " ms"); 
          INFO("Evaluating PDF         : " << t_eval         << " ms"); 
          INFO("Accept/reject          : " << t_acceptReject << " ms"); 
          INFO("Efficiency             = " << double(N) * 100. / totalGenerated   << " %");
        }
        template <typename pdf_t, typename = typename std::enable_if<!std::is_integral<pdf_t>::value>::type>
        EventList generate(pdf_t& pdf, const size_t& nEvents )
          {
            eventlist_t evts( m_eventType );
            fillEventList( pdf, evts, nEvents );
            EventList output( m_eventType );
            for( const auto& event : evts ) output.emplace_back( event );
            return output; 
          }
        EventList generate(const size_t& nEvents)
        {
          eventlist_t evts( m_eventType );
          evts.resize( nEvents ); 
          fillEventListPhaseSpace( evts, nEvents);
          EventList output( m_eventType );
          for( const auto& event : evts ) output.emplace_back( event );
          return output; 
        }


        //************SHENGHUI ADDED************//
	      template <typename eventList_t, typename pdf_t> void fillEventsUsingLambda( pdf_t& pdf, eventList_t& list, const size_t& N)
      	{
		      fillEventsUsingLambda( pdf, list, N, nullptr);
	      }
	      template <typename eventList_t, typename pdf_t, typename cut_t> void fillEventsUsingLambda( pdf_t& pdf, eventList_t& list, const size_t& N, cut_t cut )
        {
          if ( m_rnd == nullptr ) 
          {
            ERROR( "Random generator not set!" );
            return;
          }
          double maxProb   = m_normalise ? 0 : 1;
          auto size0       = list.size();
          double totalGenerated = 0; 
          // pdf.reset( true );
          ProgressBar pb(60, detail::trimmedString(__PRETTY_FUNCTION__) );
          ProfileClock t_phsp, t_eval, t_acceptReject, t_total;
          std::vector<bool> efficiencyReport(m_generatorBlock,false); 
          while ( list.size() - size0 < N ) {
            t_phsp.start();
            eventlist_t mc( m_eventType );
            mc.resize(m_generatorBlock);
            fillEventListPhaseSpace(mc, m_generatorBlock);
            t_phsp.stop();
            t_eval.start();
            // pdf.setEvents( mc );
            // pdf.prepare();
            auto previousSize = list.size();
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for ( size_t block=0; block < mc.nBlocks(); ++block )
            { 
              mc.setGenPDF(block, pdf(mc[block]) / mc.genPDF(block) );
            }
            maxProb = maxProb == 0 ? 1.2 * getMax(mc, pdf) : maxProb; 

            t_eval.stop();
            t_acceptReject.start(); 
            totalGenerated += mc.size();
            for(const auto& event : mc)
            { 
              if ( event.genPdf()  > maxProb ) {
                std::cout << std::endl; 
                WARNING( "PDF value exceeds norm value: " << event.genPdf() << " > " << maxProb );
              }
              if ( event.genPdf() > maxProb * m_rnd->Rndm() ){
                list.push_back(event);
                list.rbegin()->setGenPdf( pdf(event) );
                efficiencyReport[event.index()] = true; 
              }
              else efficiencyReport[event.index()] = false; 
              if ( list.size() - size0 == N ) break; 
            }
            t_acceptReject.stop(); 

            // m_gps.provideEfficiencyReport( efficiencyReport );
            double efficiency = 100. * ( list.size() - previousSize ) / (double)m_generatorBlock;
            pb.print( double(list.size()) / double(N), " ε[gen] = " + mysprintf("%.4f",efficiency) + "% , " + std::to_string(int(t_total.count()/1000.))  + " seconds" );
            if ( list.size() == previousSize ) {
              ERROR( "No events generated, PDF: is likely to be malformed" );
              break;
            }
          } 
          pb.finish();
          t_total.stop();
          INFO("Generated " << N << " events in " << t_total << " ms");
          INFO("Generating phase space : " << t_phsp         << " ms"); 
          INFO("Evaluating PDF         : " << t_eval         << " ms"); 
          INFO("Accept/reject          : " << t_acceptReject << " ms"); 
          INFO("Efficiency             = " << double(N) * 100. / totalGenerated   << " %");
        }

        // Fill two event lists at once
        template <typename eventList_t, typename pdf_t> void fill2EventsUsingLambda( pdf_t& pdf, eventList_t& list1, eventList_t& list2, const size_t& N)
				{
					fill2EventsUsingLambda( pdf, list1, list2, N, nullptr);
				}
				template <typename eventList_t, typename pdf_t, typename cut_t> void fill2EventsUsingLambda( pdf_t& pdf, eventList_t& list1, eventList_t& list2, const size_t& N, cut_t cut )
        {
          if ( m_rnd == nullptr ) 
          {
            ERROR( "Random generator not set!" );
            return;
          }
          double maxProb   = m_normalise ? 0 : 1;
          auto size0       = list1.size();
          double totalGenerated = 0; 
          // pdf.reset( true );
          ProgressBar pb(60, detail::trimmedString(__PRETTY_FUNCTION__) );
          ProfileClock t_phsp, t_eval, t_acceptReject, t_total;
          std::vector<bool> efficiencyReport(m_generatorBlock,false); 

          while ( list1.size() - size0 < N ) {
            t_phsp.start();
            eventlist_t mc1( list1.eventType() );
            eventlist_t mc2( list2.eventType() );
            mc1.resize(m_generatorBlock);
            mc2.resize(m_generatorBlock);
            fillEventListPhaseSpace(mc1, m_generatorBlock);
            fillEventListPhaseSpace(mc2, m_generatorBlock); // problem here?
            t_phsp.stop();
            t_eval.start();
            // pdf.setEvents( mc );
            // pdf.prepare();
            auto previousSize = list1.size();
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for ( size_t block=0; block < mc1.nBlocks(); ++block )
            { 
              mc1.setGenPDF(block, pdf(mc1[block], mc2[block]) / mc1.genPDF(block) );
              mc2.setGenPDF(block, pdf(mc1[block], mc2[block]) / mc2.genPDF(block) );
            }
            maxProb = maxProb == 0 ? 1.2 * getMax(mc1, pdf) : maxProb; // same pdf for both mc1 and 2 so only need tod o this once

            t_eval.stop();
            t_acceptReject.start(); 
            totalGenerated += mc1.size();
            // for(const auto& event : mc)
            for(int i=0; i<mc1.size(); i++)
            { 
              Event event1 = mc1[i];
              Event event2 = mc2[i];
              if ( event1.genPdf()  > maxProb || event2.genPdf() > maxProb) {
                std::cout << std::endl; 
                WARNING( "PDF value exceeds norm value: " << event1.genPdf() << " or " << event2.genPdf() << " > " << maxProb );
              }
              if ( event1.genPdf() > maxProb * m_rnd->Rndm() ){ // only need event 1 as there is only one joint pdf for both events
                list1.push_back(event1);
                list2.push_back(event2);
                list1.rbegin()->setGenPdf( pdf(event1, event2) );
                // list2.rbegin()->setGenPdf( pdf(event1, event2) );
                efficiencyReport[i] = true; 
              }
              else efficiencyReport[i] = false; 
              if ( list1.size() - size0 == N ) break; 
            }
            t_acceptReject.stop(); 

            // m_gps.provideEfficiencyReport( efficiencyReport );
            double efficiency = 100. * ( list1.size() - previousSize ) / (double)m_generatorBlock;
            pb.print( double(list1.size()) / double(N), " ε[gen] = " + mysprintf("%.4f",efficiency) + "% , " + std::to_string(int(t_total.count()/1000.))  + " seconds" );
            if ( list1.size() == previousSize ) {
              ERROR( "No events generated, PDF: is likely to be malformed" );
              break;
            }
          } 
          pb.finish();
          t_total.stop();
          INFO("Generated " << N << " events in " << t_total << " ms");
          INFO("Generating phase space : " << t_phsp         << " ms"); 
          INFO("Evaluating PDF         : " << t_eval         << " ms"); 
          INFO("Accept/reject          : " << t_acceptReject << " ms"); 
          INFO("Efficiency             = " << double(N) * 100. / totalGenerated   << " %");
        }
    };



  template <class FCN> class PDFWrapper 
  {
    public:
      void prepare(){};
      void setEvents( AmpGen::EventList& /*evts*/ ){};
      double prob_unnormalised( const AmpGen::Event& evt ) const { return m_fcn(evt); }
      explicit PDFWrapper( const FCN& fcn ) : m_fcn(fcn) {}
      size_t size() const { return 0; }
      void reset( const bool& /*flag*/ = false ){};

    private:
      FCN m_fcn;
  };

  template <class FCN> PDFWrapper<FCN> make_pdf(const FCN& fcn){ return PDFWrapper<FCN>(fcn) ; }

  /** @function PyGenerate

    Wrapper around the a phase space generator from a stringy event type to be used with python / numpy.
    */
  extern "C" void PyGenerate( const char* eventType, double* out, const unsigned int size );
} // namespace AmpGen
#endif
