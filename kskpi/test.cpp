template <typename PDF>
FitResult* doFit( PDF&& pdf, EventList_type& data, EventList_type& mc, MinuitParameterSet& MPS )
{
  auto time_wall = std::chrono::high_resolution_clock::now();
  auto time      = std::clock();

  pdf.setEvents( data );
  auto& signalPdf = std::get<0>(pdf.pdfs());
  Minimiser mini( signalPdf, &MPS );
  mini.doFit();
  FitResult* fr = new FitResult(mini);

  auto twall_end  = std::chrono::high_resolution_clock::now();
  double time_cpu = ( std::clock() - time ) / (double)CLOCKS_PER_SEC;
  double tWall    = std::chrono::duration<double, std::milli>( twall_end - time_wall ).count();
  INFO( "Wall time = " << tWall / 1000. );
  INFO( "CPU  time = " << time_cpu );
  fr->print();
  return fr;
}
auto fr = doFit(make_pdf<EventList_type>(sig,bkg), events, eventsMC, MPS)

