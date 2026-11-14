// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SymmetricTetrahedronRule.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::quadrature {

namespace {

/**
 * Orbits of a point under the permutations of the vertices, as patterns of barycentric
 * coordinates. With the parameters a, b, c:
 *   S1   AAAA  A = 1/4
 *   S4   AAAB  A = a, B = 1 - 3a
 *   S6   AABB  A = a, B = 1/2 - a                            (a < 1/4)
 *   S12  AABC  A = a, B = b, C = 1 - 2a - b                  (b < C)
 *   S24  ABCD  A = a, B = b, C = c, D = 1 - a - b - c        (a < b < c < D)
 */
enum class OrbitKind { S1, S4, S6, S12, S24 };

struct Orbit {
  OrbitKind kind;
  //! The weight of each point of the orbit, relative to the volume.
  double weight;
  std::array<double, 3> parameters;
};

struct Rule {
  std::size_t degree;
  std::vector<Orbit> orbits;
};

/**
 * The rules of quad_rules.SymmetricTet in SeisSol/matrices up to degree 14, as orbits sorted by
 * kind and then by parameters; values rounded to 20 significant digits. Sources (files in
 * src/quad_rules/tetrahedron there):
 *   degrees 1, 2, 3, 5, 6, 7, 8: F. D. Witherden, P. E. Vincent, On the identification of
 *     symmetric quadrature rules for finite element methods, Comput. Math. Appl. 69(10),
 *     1232-1241, 2015, via PyFR (BSD 3-Clause License, Copyright (c) 2012-2026 Imperial College
 *     London)
 *   degree 10: computed for SeisSol by a search over orbit structures
 *   degree 11: J. Jaskowiec, N. Sukumar, High-order symmetric cubature rules for tetrahedra and
 *     pyramids, Int. J. Numer. Methods Eng. 122(1), 148-171, 2021 (CC BY 4.0)
 *   degree 13: F. D. Witherden, via PyFR (as above)
 *   degrees 12, 14: computed for SeisSol from the rules of Jaskowiec and Sukumar of degree 13
 *     and 14 by removing and merging orbits
 */
const std::vector<Rule>& rules() {
  // one orbit per line; designators would spread every orbit over three lines
  // NOLINTBEGIN(modernize-use-designated-initializers)
  static const std::vector<Rule> Rules = {
      // witherden-vincent/witherden-vincent-n1-d1.txt
      {1,
       {
           {OrbitKind::S1, 1.0, {}},
       }},
      // witherden-vincent/witherden-vincent-n4-d2.txt
      {2,
       {
           {OrbitKind::S4, 0.25, {0.13819660112501051518}},
       }},
      // witherden-vincent/witherden-vincent-n8-d3.txt
      {3,
       {
           {OrbitKind::S4, 0.11378215746291264293, {0.10804724989842860412}},
           {OrbitKind::S4, 0.13621784253708735707, {0.32816330251638168679}},
       }},
      // witherden-vincent/witherden-vincent-n14-d5.txt
      {5,
       {
           {OrbitKind::S4, 0.073493043116361949544, {0.092735250310891226402}},
           {OrbitKind::S4, 0.1126879257180158508, {0.3108859192633006098}},
           {OrbitKind::S6, 0.042546020777081466438, {0.045503704125649649492}},
       }},
      // witherden-vincent/witherden-vincent-n24-d6.txt
      {6,
       {
           {OrbitKind::S4, 0.010077211055320642948, {0.040673958534611353116}},
           {OrbitKind::S4, 0.0399227502581674921, {0.21460287125915202929}},
           {OrbitKind::S4, 0.055357181543654722095, {0.32233789014227551034}},
           {OrbitKind::S12,
            0.048214285714285714286,
            {0.063661001875017525299, 0.26967233145831580803}},
       }},
      // witherden-vincent/witherden-vincent-n35-d7.txt
      {7,
       {
           {OrbitKind::S1, 0.095485289464130848861, {}},
           {OrbitKind::S4, 0.042329581209967029076, {0.31570114977820279942}},
           {OrbitKind::S6, 0.031896927832857579934, {0.050489822598396368763}},
           {OrbitKind::S12,
            0.0081107708299033415661,
            {0.021265472541483245989, 0.1466388138184849469}},
           {OrbitKind::S12,
            0.03720713072833462137,
            {0.18883383102600104774, 0.047160700360997881044}},
       }},
      // witherden-vincent/witherden-vincent-n46-d8.txt
      {8,
       {
           {OrbitKind::S4, 0.0075252561535401997807, {0.042316543684767281624}},
           {OrbitKind::S4, 0.026426650908408831404, {0.10795272496221084848}},
           {OrbitKind::S4, 0.052031747563738532715, {0.18510948778258657028}},
           {OrbitKind::S4, 0.041763782856934897707, {0.31418170912403899536}},
           {OrbitKind::S6, 0.036280930261308823001, {0.064408671416169794772}},
           {OrbitKind::S12,
            0.0071569028908444331502,
            {0.021433930127130574844, 0.23966807631943052702}},
           {OrbitKind::S12,
            0.015453486150960334814,
            {0.20413933387602912043, 0.0079239539457973532127}},
       }},
      // orbit-search/orbit-search-n79-d10.txt
      {10,
       {
           {OrbitKind::S1, 0.045744342446789630054, {}},
           {OrbitKind::S4, 0.00055346716868564207782, {0.010638211251437933813}},
           {OrbitKind::S4, 0.01092468847059769422, {0.11422900089056276921}},
           {OrbitKind::S4, 0.025693972060698947718, {0.31273589377082560697}},
           {OrbitKind::S6, 0.00055953733740315898993, {0.016432313689786628482}},
           {OrbitKind::S12,
            0.0061585332282130537542,
            {0.030018728177575208716, 0.12782018682468132107}},
           {OrbitKind::S12,
            0.010445427028362449501,
            {0.034305206188598863068, 0.33305715558157544454}},
           {OrbitKind::S12,
            0.025140755309605840298,
            {0.12345767363800584155, 0.28185872506094552653}},
           {OrbitKind::S12,
            0.01331824037907913148,
            {0.17396983292677239303, 0.022886836789337276292}},
           {OrbitKind::S12,
            0.011787870948811381628,
            {0.40991377504957833642, 0.014698521078267115984}},
       }},
      // jaskowiec-sukumar-symmetric/jaskowiec-sukumar-n110-d11.txt
      {11,
       {
           {OrbitKind::S4, 0.0033808695306217, {0.0326328439651818}},
           {OrbitKind::S4, 0.0146893435861307, {0.2993130923993425}},
           {OrbitKind::S6, 0.0016532536103758, {0.0057019051227004}},
           {OrbitKind::S6, 0.0193236403477585, {0.0987846241416251}},
           {OrbitKind::S6, 0.0167667686772256, {0.1742044154846892}},
           {OrbitKind::S12, 0.0020362784913048, {0.0130038704855902, 0.1702280734892642}},
           {OrbitKind::S12, 0.0082076137529167, {0.0429354011719257, 0.2958693878670328}},
           {OrbitKind::S12, 0.0097999712470329, {0.1231206549747224, 0.0264283437041717}},
           {OrbitKind::S12, 0.0151049873666635, {0.1844240320470778, 0.0912518190501428}},
           {OrbitKind::S12, 0.0127423926661064, {0.2591694074642826, 0.0262258654494753}},
           {OrbitKind::S24,
            0.0052734270596891,
            {0.0109775723659603, 0.1124285605067612, 0.341451866195819}},
       }},
      // reduced-symmetric/reduced-symmetric-n128-d12.txt
      {12,
       {
           {OrbitKind::S4, 0.0012213236394122397463, {0.023245033439622245205}},
           {OrbitKind::S4, 0.0041683480772906546102, {0.067818444853287975615}},
           {OrbitKind::S4, 0.01646452592137550917, {0.19227542350347669917}},
           {OrbitKind::S4, 0.01163914863855245242, {0.2995946664293765296}},
           {OrbitKind::S4, 0.0095203374361476411075, {0.32856535244020992646}},
           {OrbitKind::S6, 0.0051896118071942558424, {0.024960354798928596994}},
           {OrbitKind::S6, 0.013621443987026121145, {0.1658414213906112664}},
           {OrbitKind::S12,
            0.0018295264187962173646,
            {0.01732335619801967238, 0.12090990286867206601}},
           {OrbitKind::S12,
            0.0048211923825337232086,
            {0.026469983069155665011, 0.27337652833520918604}},
           {OrbitKind::S12,
            0.013252225334128021319,
            {0.088711963156791796884, 0.21548879191713861988}},
           {OrbitKind::S12,
            0.0053819114555618086685,
            {0.12214832216526455215, 0.01617980740530676795}},
           {OrbitKind::S12,
            0.013526338349954384571,
            {0.2295849225845133027, 0.055181875740560524005}},
           {OrbitKind::S12,
            0.011903953152938325467,
            {0.41008417834480802056, 0.053669171802587585089}},
           {OrbitKind::S24,
            0.0044373818856922492783,
            {0.0047101024834800475894, 0.13345554883393696269, 0.30869155100705005208}},
       }},
      // witherden/witherden-n145-d13.txt
      {13,
       {
           {OrbitKind::S1, 0.027909817545077277362, {}},
           {OrbitKind::S4, 0.00066797126602980607763, {0.01830475748619279897}},
           {OrbitKind::S4, 0.019307298608679023068, {0.17901508263002278126}},
           {OrbitKind::S4, 0.0066968089575574080063, {0.32961585375444814701}},
           {OrbitKind::S6, 0.0025366658199614153096, {0.01574108080395255591}},
           {OrbitKind::S6, 0.0085645146236781589377, {0.062306200622718381386}},
           {OrbitKind::S12,
            0.0019633582546351872041,
            {0.01685802508195034297, 0.24966122768071147328}},
           {OrbitKind::S12,
            0.0025267914252938963788,
            {0.023525718244859602976, 0.10234763590487778111}},
           {OrbitKind::S12,
            0.0059936450115461083784,
            {0.069441895049546465147, 0.15254189508030244848}},
           {OrbitKind::S12,
            0.011317403885189328579,
            {0.093700527282147662399, 0.24914221761321528683}},
           {OrbitKind::S12,
            0.0032210831532635021363,
            {0.125885360164042453, 0.0081232940042499605902}},
           {OrbitKind::S12,
            0.010795886453342852088,
            {0.21795527054773767475, 0.028798833628512335582}},
           {OrbitKind::S12,
            0.017511763976330354249,
            {0.35466345547278388496, 0.090658178843814455091}},
           {OrbitKind::S12,
            0.0057765015420767012478,
            {0.41668928765703838516, 0.014566569709145578335}},
           {OrbitKind::S24,
            0.0037298991684952152258,
            {0.011242443361469569299, 0.077946562231807870883, 0.30415843367637258712}},
       }},
      // reduced-symmetric/reduced-symmetric-n188-d14.txt
      {14,
       {
           {OrbitKind::S4, 0.00012669861498529424145, {0.0071271917846273425924}},
           {OrbitKind::S4, 0.0031556430453557847132, {0.060181316144895010593}},
           {OrbitKind::S4, 0.0073807550439583712158, {0.10506704683668688399}},
           {OrbitKind::S4, 0.012850438163312283702, {0.17094032168444347954}},
           {OrbitKind::S4, 0.0035434886730392892481, {0.33001686227002727498}},
           {OrbitKind::S6, 0.00057376145133211772357, {0.0047726292184414800979}},
           {OrbitKind::S6, 0.0029173681344473621933, {0.031906951512464572213}},
           {OrbitKind::S6, 0.013404394023084167823, {0.10043778348655375647}},
           {OrbitKind::S6, 0.014306359644459105898, {0.18273833873917376154}},
           {OrbitKind::S12,
            0.00087810368182705930623,
            {0.012623387734500326889, 0.12738769333130849825}},
           {OrbitKind::S12,
            0.0019345417048007533845,
            {0.016653990541828016459, 0.28411430763669068822}},
           {OrbitKind::S12,
            0.00083267169654926401232,
            {0.04931011942255638003, 0.006027821334791865917}},
           {OrbitKind::S12,
            0.004318631869500129221,
            {0.21375251944142869638, 0.0092879473773393915062}},
           {OrbitKind::S12,
            0.013122721976718622152,
            {0.25122195380394712836, 0.07775995191888354886}},
           {OrbitKind::S12,
            0.0070822514075175054762,
            {0.39445742405546280419, 0.020385785977262600074}},
           {OrbitKind::S24,
            0.0041247421820764962157,
            {0.01508741786999424775, 0.082589174362511825503, 0.34924776400456941214}},
           {OrbitKind::S24,
            0.0037735955723057893958,
            {0.017923115648495079027, 0.078836945793028871424, 0.16954602309105399161}},
           {OrbitKind::S24,
            0.007373893007055188683,
            {0.060758751957025510785, 0.12274782493550955677, 0.24114955109756799626}},
       }},
  };
  // NOLINTEND(modernize-use-designated-initializers)
  return Rules;
}

//! The barycentric coordinates of the points of an orbit: the distinct permutations of its
//! pattern, in lexicographic order of the letters.
std::vector<std::array<double, 4>> orbitPoints(const Orbit& orbit) {
  const auto& p = orbit.parameters;
  std::array<int, 4> letters{};
  std::array<double, 4> values{};
  switch (orbit.kind) {
  case OrbitKind::S1:
    letters = {0, 0, 0, 0};
    values = {0.25, 0, 0, 0};
    break;
  case OrbitKind::S4:
    letters = {0, 0, 0, 1};
    values = {p[0], 1 - 3 * p[0], 0, 0};
    break;
  case OrbitKind::S6:
    letters = {0, 0, 1, 1};
    values = {p[0], 0.5 - p[0], 0, 0};
    break;
  case OrbitKind::S12:
    letters = {0, 0, 1, 2};
    values = {p[0], p[1], 1 - 2 * p[0] - p[1], 0};
    break;
  case OrbitKind::S24:
    letters = {0, 1, 2, 3};
    values = {p[0], p[1], p[2], 1 - p[0] - p[1] - p[2]};
    break;
  }
  std::vector<std::array<double, 4>> points;
  do {
    points.push_back(
        {values[letters[0]], values[letters[1]], values[letters[2]], values[letters[3]]});
  } while (std::ranges::next_permutation(letters).found);
  return points;
}

} // namespace

std::pair<std::vector<std::array<double, 3>>, std::vector<double>>
    symmetricTetrahedronRule(std::size_t degree) {
  const auto& all = rules();
  const auto rule =
      std::ranges::find_if(all, [&](const Rule& candidate) { return candidate.degree >= degree; });
  if (rule == all.end()) {
    logError() << "No symmetric tetrahedral quadrature rule of degree" << degree
               << "is available; the highest is" << MaxSymmetricTetrahedronDegree;
  }

  std::pair<std::vector<std::array<double, 3>>, std::vector<double>> result;
  for (const auto& orbit : rule->orbits) {
    for (const auto& lambda : orbitPoints(orbit)) {
      // the vertex (0,0,0) carries lambda[0]
      result.first.push_back({lambda[1], lambda[2], lambda[3]});
      result.second.push_back(orbit.weight / 6.0);
    }
  }
  return result;
}

} // namespace seissol::quadrature
