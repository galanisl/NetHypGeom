/*
 Copyright (c) 2012 Thijs (M.A.) van den Berg, http://www.sitmo.com/
 (See https://code.google.com/p/stdfin/source/browse/trunk/LICENSE.txt)
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of
 this software and associated documentation files (the "Software"), to deal in
 the Software without restriction, including without limitation the rights to
 use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 of the Software, and to permit persons to whom the Software is furnished to do
 so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#ifndef STDFIN_NUMERICAL_INTEGRATION_GAUSS_LEGENDRE_HPP
#define STDFIN_NUMERICAL_INTEGRATION_GAUSS_LEGENDRE_HPP

namespace stdfin {

namespace ni {
//---------------------------------------------------------------------
// 8 point rule
//---------------------------------------------------------------------
const double gauss_legendre_w8[] = {
  0.362683783378362,  0.313706645877887, 0.222381034453374, 0.101228536290376};
const double gauss_legendre_x8[] = {
  0.183434642495650, 0.525532409916329, 0.796666477413627, 0.960289856497536};

//---------------------------------------------------------------------
// 16 point rule
//---------------------------------------------------------------------
const double gauss_legendre_w16[] = {
  0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
  0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541};
const double gauss_legendre_x16[] = {
  0.0950125098376374, 0.2816035507792589, 0.4580167776572274, 0.6178762444026437,
  0.7554044083550030, 0.8656312023878317, 0.9445750230732326, 0.9894009349916499};

//---------------------------------------------------------------------
// 24 point rule
//---------------------------------------------------------------------
const double gauss_legendre_w24[] = {
  0.1279381953467522, 0.1258374563468283, 0.1216704729278034, 0.1155056680537256,
  0.1074442701159656, 0.0976186521041139, 0.0861901615319533, 0.0733464814110803,
  0.0592985849154368, 0.0442774388174198, 0.0285313886289337, 0.0123412297999872};
const double gauss_legendre_x24[] = {
  0.0640568928626056, 0.1911188674736163, 0.3150426796961634, 0.4337935076260451,
  0.5454214713888395, 0.6480936519369756, 0.7401241915785544, 0.8200019859739029,
  0.8864155270044010, 0.9382745520027328, 0.9747285559713095, 0.9951872199970214};

//---------------------------------------------------------------------
// 32 point rule
//---------------------------------------------------------------------
const double gauss_legendre_w32[] = {
  0.0965400885147278, 0.0956387200792749, 0.0938443990808046, 0.0911738786957639,
  0.0876520930044038, 0.0833119242269467, 0.0781938957870703, 0.0723457941088485,
  0.0658222227763618, 0.0586840934785355, 0.0509980592623762, 0.0428358980222267,
  0.0342738629130214, 0.0253920653092621, 0.0162743947309057, 0.0070186100094701};
const double gauss_legendre_x32[] = {
  0.0483076656877383, 0.1444719615827960, 0.2392873622521370, 0.3318686022821270,
  0.4213512761306350, 0.5068999089322290, 0.5877157572407620, 0.6630442669302150,
  0.7321821187402890, 0.7944837959679420, 0.8493676137325700, 0.8963211557660520,
  0.9349060759377390, 0.9647622555875060, 0.9856115115452680, 0.9972638618494810};

//---------------------------------------------------------------------
// 48 point rule
//---------------------------------------------------------------------
const double gauss_legendre_w48[] = {
  0.064737696812683923, 0.064466164435950082, 0.063924238584648187, 0.063114192286254026,
  0.062039423159892664, 0.060704439165893880, 0.059114839698395636, 0.057277292100403216,
  0.055199503699984163, 0.052890189485193667, 0.050359035553854475, 0.047616658492490475,
  0.044674560856694280, 0.041545082943464749, 0.038241351065830706, 0.034777222564770439,
  0.031167227832798089, 0.027426509708356948, 0.023570760839324379, 0.019616160457355528,
  0.015579315722943849, 0.011477234579234539, 0.007327553901276262, 0.003153346052305839};
const double gauss_legendre_x48[] = {
  0.032380170962869362, 0.097004699209462699, 0.161222356068891718, 0.224763790394689061,
  0.287362487355455577, 0.348755886292160738, 0.408686481990716730, 0.466902904750958405,
  0.523160974722233034, 0.577224726083972704, 0.628867396776513624, 0.677872379632663905,
  0.724034130923814655, 0.767159032515740339, 0.807066204029442627, 0.843588261624393531,
  0.876572020274247886, 0.905879136715569673, 0.931386690706554333, 0.952987703160430861,
  0.970591592546247250, 0.984124583722826858, 0.993530172266350758, 0.998771007252426119};
}


#define INTEGRATE_GAUSS_LEGENDRE(N,HN)                                     \
template <class F>                                                         \
double gauss_legendre_##N(F& f, double a, double b) {                      \
  double sum = 0;                                                          \
  double rescale = 0.5 * (b - a);                                          \
  double center  = 0.5 * (b + a);                                          \
  for (int i = 0; i < HN; ++i) {                                           \
    double dx = rescale * ni::gauss_legendre_x##N[i];                      \
    sum += ni::gauss_legendre_w##N[i]*( f(center + dx) + f(center - dx) ); \
  }                                                                        \
  return sum * rescale;                                                    \
}

INTEGRATE_GAUSS_LEGENDRE(8,4)
  INTEGRATE_GAUSS_LEGENDRE(16,8)
  INTEGRATE_GAUSS_LEGENDRE(24,12)
  INTEGRATE_GAUSS_LEGENDRE(32,16)
  INTEGRATE_GAUSS_LEGENDRE(48,24)
  
} // namespace 
#endif

