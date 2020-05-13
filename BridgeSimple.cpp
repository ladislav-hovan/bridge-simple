/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"
#include "colvar/Colvar.h"
#include <iostream>

namespace PLMD {
namespace colvar {

	class BridgeSimple : public Colvar {
	private:
		std::vector<AtomNumber> ga_lista, gb_lista, b_lista, full_lista, red_lista;
		unsigned b_start, b_finish;
		SwitchingFunction sf1;
		SwitchingFunction sf2;
		bool doneigh = false, firsttime = true, invalidateList = true;
		double nl_cut = 0.0;
		int nl_st = 0;
	public:
		static void registerKeywords(Keywords& keys);
		explicit BridgeSimple(const ActionOptions&);
		// Active methods:
		//double compute(unsigned j, unsigned k, std::vector<Vector>& deriv, Tensor& virial) const;
		virtual void prepare();
		virtual void calculate();
		bool isPeriodic() { return false; }
	};

	PLUMED_REGISTER_ACTION(BridgeSimple, "BRIDGE_SIMPLE")

	void BridgeSimple::registerKeywords(Keywords& keys) {
		Colvar::registerKeywords(keys);
		keys.add("atoms-2", "BRIDGING_ATOMS", "The list of atoms that can form the bridge between the two interesting parts "
			"of the structure.");
		keys.add("atoms-2", "GROUPA", "The list of atoms that are in the first interesting part of the structure");
		keys.add("atoms-2", "GROUPB", "The list of atoms that are in the second interesting part of the structure");
		keys.add("optional", "SWITCH", "The parameters of the two \\ref switchingfunction in the above formula");
		keys.add("optional", "SWITCHA", "The \\ref switchingfunction on the distance between bridging atoms and the atoms in "
			"group A");
		keys.add("optional", "SWITCHB", "The \\ref switchingfunction on the distance between the bridging atoms and the atoms in "
			"group B");
		keys.addFlag("NLIST", false, "Use a neighbor list to speed up the calculation");
		keys.add("optional", "NL_CUTOFF", "The cutoff for the neighbor list");
		keys.add("optional", "NL_STRIDE", "The frequency with which we are updating the atoms in the neighbor list");
	}

	BridgeSimple::BridgeSimple(const ActionOptions&ao) :
		PLUMED_COLVAR_INIT(ao)
	{
		addValueWithDerivatives(); setNotPeriodic();

		// Read in the atoms
		// std::vector<AtomNumber> ga_lista, gb_lista, b_lista;
		parseAtomList("GROUPA", ga_lista);
		parseAtomList("GROUPB", gb_lista);
		parseAtomList("BRIDGING_ATOMS", b_lista);

		// Neighbour list stuff
		parseFlag("NLIST", doneigh);
		if (doneigh) {
			parse("NL_CUTOFF", nl_cut);
			if (nl_cut <= 0.0) error("NL_CUTOFF should be explicitly specified and positive");
			parse("NL_STRIDE", nl_st);
			if (nl_st <= 0) error("NL_STRIDE should be explicitly specified and positive");
		}

		// Read the switches
		std::string sfinput, errors; parse("SWITCH", sfinput);
		if (sfinput.length() > 0) {
			sf1.set(sfinput, errors);
			if (errors.length() != 0) error("problem reading SWITCH keyword : " + errors);
			sf2.set(sfinput, errors);
			if (errors.length() != 0) error("problem reading SWITCH keyword : " + errors);
		}
		else {
			parse("SWITCHA", sfinput);
			if (sfinput.length() > 0) {
				sf1.set(sfinput, errors);
				if (errors.length() != 0) error("problem reading SWITCHA keyword : " + errors);
				sfinput.clear(); parse("SWITCHB", sfinput);
				if (sfinput.length() == 0) error("found SWITCHA keyword without SWITCHB");
				sf2.set(sfinput, errors);
				if (errors.length() != 0) error("problem reading SWITCHB keyword : " + errors);
			}
			else {
				error("missing definition of switching functions");
			}
		}
		log.printf("  distance between bridging atoms and atoms in GROUPA must be less than %s\n", sf1.description().c_str());
		log.printf("  distance between bridging atoms and atoms in GROUPB must be less than %s\n", sf2.description().c_str());

		// Request all the atoms
		full_lista = ga_lista;
		full_lista.insert(full_lista.end(), gb_lista.begin(), gb_lista.end());
		b_start = full_lista.size();
		full_lista.insert(full_lista.end(), b_lista.begin(), b_lista.end());
		b_finish = full_lista.size();
		red_lista = full_lista;
		requestAtoms(full_lista);

		// And check everything has been read in correctly
		checkRead();
	}

	void BridgeSimple::prepare() {
		if (nl_st > 0) {
			if (firsttime == true || getStep() % nl_st == 0) {
				requestAtoms(full_lista);
				b_finish = full_lista.size();
				invalidateList = true;
				firsttime = false;
			}
			else {
				requestAtoms(red_lista);
				b_finish = red_lista.size();
				invalidateList = false;
				if (getExchangeStep()) error("Neighbor lists should be updated on exchange steps - choose a NL_STRIDE which divides the exchange stride!");
			}
			if (getExchangeStep()) firsttime = true;
		}
	}

	//double BridgeSimple::compute(unsigned j, unsigned k, std::vector<Vector>& deriv, Tensor& virial) const {
	//	double tot = 0.0;
	//	for (unsigned i = b_start; i < b_finish; ++i) {
	//		Vector dij = pbcDistance(getPosition(i), getPosition(j));
	//		double dw1, w1 = sf1.calculateSqr(dij.modulo2(), dw1);
	//		if (w1 == 0.0 && dw1 == 0.0)
	//			continue;
	//		Vector dik = pbcDistance(getPosition(i), getPosition(k));
	//		double dw2, w2 = sf2.calculateSqr(dik.modulo2(), dw2);
	//		if (w2 == 0.0 && dw2 == 0.0)
	//			continue;

	//		tot += w1 * w2;
	//		// And finish the calculation
	//		deriv[j] += w2 * dw1 * dij;
	//		deriv[k] += w1 * dw2 * dik;
	//		deriv[i] -= (w1 * dw2 * dik + w2 * dw1 * dij);
	//		virial += (w1 * (-dw2)*Tensor(dik, dik) + w2 * (-dw1)*Tensor(dij, dij));
	//	}
	//	return tot;
	//}

	void BridgeSimple::calculate() {
		std::vector<Vector> deriv(getNumberOfAtoms());

		if (nl_st > 0 && invalidateList) {
			red_lista.erase(red_lista.begin() + b_start, red_lista.end());
			double limit = nl_cut * nl_cut;
			for (unsigned i = b_start; i < b_finish; ++i)
			{
				bool close = false;
				for (unsigned j = 0; j < ga_lista.size(); ++j)
					if (pbcDistance(getPosition(i), getPosition(j)).modulo2() <= limit)
					{
						close = true;
						break;
					}
				if (close)
					for (unsigned k = ga_lista.size(); k < b_start; ++k)
						if (pbcDistance(getPosition(i), getPosition(k)).modulo2() <= limit)
						{
							red_lista.push_back(b_lista[i - b_start]);
							break;
						}
			}
		}

		double value = 0.0;
		Tensor virial;

		for (unsigned i = b_start; i < b_finish; ++i) {
			for (unsigned j = 0; j < ga_lista.size(); ++j) {
				Vector dij = pbcDistance(getPosition(i), getPosition(j));
				double dw1, w1 = sf1.calculateSqr(dij.modulo2(), dw1);
				if (w1 == 0.0 && dw1 == 0.0)
					continue;

				for (unsigned k = ga_lista.size(); k < b_start; ++k) {
					Vector dik = pbcDistance(getPosition(i), getPosition(k));
					double dw2, w2 = sf2.calculateSqr(dik.modulo2(), dw2);
					if (w2 == 0.0 && dw2 == 0.0)
						continue;

					value += w1 * w2;
					// And finish the calculation
					deriv[j] += w2 * dw1 * dij;
					deriv[k] += w1 * dw2 * dik;
					deriv[i] -= (w1 * dw2 * dik + w2 * dw1 * dij);
					virial += (w1 * (-dw2)*Tensor(dik, dik) + w2 * (-dw1)*Tensor(dij, dij));
				}
			}
		}

		for (unsigned i = 0; i < deriv.size(); ++i)
		{
			setAtomsDerivatives(i, deriv[i]); 
			//std::cout << deriv[i][0] << " ";
		}
		//std::cout << std::endl;
		setValue(value);
		setBoxDerivatives(virial);
	}
}
}