/*----------------------
Copyright (C): OpenGATE Collaboration

This software is distributed under the terms
of the GNU Lesser General  Public Licence (LGPL)
See LICENSE.md for further details
----------------------*/

#include "G4EmParameters.hh"
#include "GateBioDoseActor.hh"
#include "GateImageWithStatistic.hh"
#include <CLHEP/Units/SystemOfUnits.h>

#define GATE_BUFFERSIZE

//-----------------------------------------------------------------------------
GateBioDoseActor::GateBioDoseActor(G4String name, G4int depth):
	GateVImageActor(std::move(name), depth),
	_currentEvent(0),
	_messenger(this),
	_alphaRef(-1),
	_betaRef(-1),
	_sobpWeight(0),
	_enableEdep(false),
	_enableDose(true),
	_enableBioDose(true),
	_enableAlphaMix(false),
	_enableBetaMix(false),
	_enableRBE(false),
	_enableUncertainties(true)
{
	GateDebugMessageInc("Actor", 4, "GateBioDoseActor() -- begin\n");
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void GateBioDoseActor::Construct() {
	GateDebugMessageInc("Actor", 4, "GateBioDoseActor -- Construct - begin\n");
	GateVImageActor::Construct();

	// Find G4_WATER
	G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
	// Find OtherMaterial
	//G4NistManager::Instance()->FindOrBuildMaterial(mOtherMaterial);
	G4EmParameters::Instance()->SetBuildCSDARange(true);

	// Enable callbacks BioDose
	EnableBeginOfRunAction(true);
	EnableEndOfRunAction(true);
	EnableBeginOfEventAction(true);
	EnableEndOfEventAction(true);
	EnablePreUserTrackingAction(false);
	EnablePostUserTrackingAction(false);
	EnableUserSteppingAction(true);

	// Outputs
	{
		G4String basename = removeExtension(mSaveFilename);
		G4String ext = getExtension(mSaveFilename);

		auto setupImage = [&](GateImageWithStatistic& image, std::string const& suffix = "") {
			SetOriginTransformAndFlagToImage(image);
			image.SetResolutionAndHalfSize(mResolution, mHalfSize, mPosition);
			image.Allocate();

			if(!suffix.empty()) {
				G4String filename = basename + "_" + suffix + "." + ext;
				image.SetFilename(filename);
			}
		};

		if(_enableEdep)     setupImage(_edepImage, "edep");
		if(_enableDose)     setupImage(_doseImage, "dose");
		if(_enableBioDose)  setupImage(_bioDoseImage, "biodose");
		if(_enableAlphaMix) setupImage(_alphaMixImage, "alphamix");
		if(_enableBetaMix)  setupImage(_betaMixImage, "betamix");
		if(_enableRBE)      setupImage(_rbeImage, "rbe");
		if(_enableUncertainties) {
			setupImage(_biodoseUncertaintyImage, "biodose_uncertainty");

			setupImage(_eventEdepImage);

			setupImage(_eventAlphaImage);
			setupImage(_squaredAlphaMixImage);

			setupImage(_eventSqrtBetaImage);
			setupImage(_squaredSqrtBetaMixImage);

			setupImage(_alphaMixSqrtBetaMixImage);
		}
	}

	ResetData();

	///////////////////////////////////////////////////////////////////////////////////////////
	//Just matrix information
	G4cout << "Memory space to store physical dose into " << mResolution.x() * mResolution.y() * mResolution.z() << " voxels has been allocated " << G4endl;

	// SOBP
	if(_sobpWeight == 0) { _sobpWeight = 1; }

	//Building the cell line information
	_dataBase = "data/" + _cellLine + "_" + _bioPhysicalModel + ".db";
	BuildDatabase();

	if(_alphaRef < 0 || _betaRef < 0)
		GateError("BioDoseActor " << GetName() << ": setAlphaRef and setBetaRef must be done");
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void GateBioDoseActor::BuildDatabase() {
	std::ifstream f(_dataBase);
	if(!f) GateError("BioDoseActor " << GetName() << ": unable to open file '" << _dataBase << "'");

	int nZ = 0;
	double prevKineticEnergy = 1;
	double prevAlpha = 1;
	double prevBeta =1;

	for(std::string line; std::getline(f, line); ) {
		std::istringstream iss(line);
		std::string firstCol;

		iss >> firstCol;

		if(firstCol == "Fragment") {
			if(nZ != 0) // prevKineticEnergy is the maximum kinetic energy for current nZ
				_energyMaxForZ[nZ] = prevKineticEnergy;

			iss >> nZ;
			prevKineticEnergy = 1;
			prevAlpha = 1;
			prevBeta = 1;
		} else if(nZ != 0) {
			double kineticEnergy = 0;
			double alpha = 0;
			double beta = 0;
			std::istringstream{firstCol} >> kineticEnergy;
			iss >> alpha;
			iss >> beta;

			auto alphaCoeff = Interpol(prevKineticEnergy, kineticEnergy, prevAlpha, alpha);
			auto sqrtBetaCoeff = Interpol(prevKineticEnergy, kineticEnergy, std::sqrt(prevBeta), std::sqrt(beta));

			// Saving the in the input databse
			Fragment fragment{nZ, kineticEnergy};
			_alphaBetaInterpolTable[fragment] = {alphaCoeff, sqrtBetaCoeff};

			prevKineticEnergy = kineticEnergy;
			prevAlpha = alpha;
			prevBeta = beta;
		} else {
			GateError("BioDoseActor " << GetName() << ": bad database format in '" << _dataBase << "'");
		}
	}

	if(nZ != 0) // last line read; prevKineticEnergy is the maximum kinetic energy for current nZ
		_energyMaxForZ[nZ] = prevKineticEnergy;
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
GateBioDoseActor::Coefficients GateBioDoseActor::Interpol(double x1, double x2, double y1, double y2) {
	//Function for a 1D linear interpolation. It returns a pair of a and b coefficients
	double a = (y2 - y1) / (x2 - x1);
	double b = y1 - x1 * a;
	return {a, b};
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void GateBioDoseActor::SaveData() {
	GateDebugMessageInc("Actor", 4, "GateBioDoseActor::SaveData() known ion events / total events: " << _eventWithKnownIonCount << " / " << _eventCount << "\n");

	auto const sqAlphaRef = _alphaRef * _alphaRef;

	std::ofstream f{"output/uncinfo", std::ios_base::trunc};

	for(auto const& [index, deposited]: _depositedMap) {
		// Alpha Beta mix (final)
		double alphaMix = 0;
		double sqrtBetaMix = 0;
		double betaMix = 0;

		if(deposited.energy != 0) {
			alphaMix = (deposited.alpha / deposited.energy);
			sqrtBetaMix = (deposited.sqrtBeta / deposited.energy);
			betaMix = sqrtBetaMix * sqrtBetaMix;
		}

		// Calculate biological dose and RBE
		double biodose  = 0;
		double rbe      = 0;

		auto const dose = deposited.dose;
		auto const sqDose = deposited.dose * deposited.dose;
		auto const sqrtDelta = std::sqrt(sqAlphaRef + 4 * _betaRef * (alphaMix * dose + betaMix * sqDose));

		if(deposited.dose > 0 && alphaMix != 0 && betaMix != 0) {
			biodose = (-_alphaRef + sqrtDelta) / (2 * _betaRef);
			rbe = biodose / dose;
		}

		if(_enableUncertainties) {
			if(deposited.dose > 0 && alphaMix != 0 && betaMix != 0 && sqrtDelta > 0 && _currentEvent > 0) {
				double n = deposited.n;
				n = _currentEvent;

				double sumSquaredAlphaMix = _squaredAlphaMixImage.GetValue(index);
				double sumSquaredSqrtBetaMix = _squaredSqrtBetaMixImage.GetValue(index);
				double sumAlphaMixSqrtBetaMix = _alphaMixSqrtBetaMixImage.GetValue(index);

				double pdBiodoseAlphaD = 1. / sqrtDelta;
				double pdBiodoseSqrtBetaD = (2. * sqrtBetaMix * dose) / sqrtDelta;
				double varAlphaMix = sumSquaredAlphaMix / n - alphaMix * alphaMix / (n*n);
				double varSqrtBetaMix = sumSquaredSqrtBetaMix / n - sqrtBetaMix * sqrtBetaMix / (n*n);
				double covAlphaMixSqrtBetaMix = sumAlphaMixSqrtBetaMix / n - (alphaMix * sqrtBetaMix) / (n*n);
				double pdBiodoseAlphaMix = dose * pdBiodoseAlphaD;
				double pdBiodoseSqrtBetaMix = dose * pdBiodoseSqrtBetaD;
				double uncertaintyBiodose = std::sqrt(
					pdBiodoseAlphaD * pdBiodoseAlphaD * varAlphaMix +
					pdBiodoseSqrtBetaD * pdBiodoseSqrtBetaD * varSqrtBetaMix +
					2. * covAlphaMixSqrtBetaMix * pdBiodoseAlphaMix * pdBiodoseSqrtBetaMix
				);

				_biodoseUncertaintyImage.SetValue(index, uncertaintyBiodose);

				{ // TODO remove for release
					f << "index: " << index << '\n';
					f << "n: " << deposited.n << '\n';
					f << "phydose: " << dose << '\n';
					f << "biodose: " << biodose << '\n';
					f << "alphaMix: " << alphaMix << '\n';
					f << "sqrtBetaMix: " << sqrtBetaMix << '\n';
					f << "sqrtDelta: " << sqrtDelta << '\n';
					f << "var_a: " << varAlphaMix << '\n';
					f << "var_b: " << varSqrtBetaMix << '\n';
					f << "pd_a: " << pdBiodoseAlphaD << '\n';
					f << "pd_b: " << pdBiodoseSqrtBetaD << '\n';
					f << "cov: " << covAlphaMixSqrtBetaMix << '\n';
					f << "alpha part:  " << pdBiodoseAlphaD * pdBiodoseAlphaD * varAlphaMix << '\n';
					f << "beta part:   " << pdBiodoseSqrtBetaD * pdBiodoseSqrtBetaD * varSqrtBetaMix << '\n';
					f << "cov part:    " << 2. * covAlphaMixSqrtBetaMix * pdBiodoseAlphaMix * pdBiodoseSqrtBetaMix << '\n';
					f << "uncertainty: " << uncertaintyBiodose << '\n';
					f << "---------------------\n";
				}
			} else {
				_biodoseUncertaintyImage.SetValue(index, 1);
			}
		}

		// Write data
		if(_enableEdep)     _edepImage.SetValue(index, deposited.energy);
		if(_enableDose)     _doseImage.SetValue(index, deposited.dose);
		if(_enableBioDose)  _bioDoseImage.SetValue(index, biodose);
		if(_enableAlphaMix) _alphaMixImage.SetValue(index, alphaMix);
		if(_enableBetaMix)  _betaMixImage.SetValue(index, betaMix);
		if(_enableRBE)      _rbeImage.SetValue(index, rbe);
	}

	GateVActor::SaveData();

	if(_enableEdep)           _edepImage.SaveData(_currentEvent);
	if(_enableDose)           _doseImage.SaveData(_currentEvent);
	if(_enableBioDose)        _bioDoseImage.SaveData(_currentEvent);
	if(_enableAlphaMix)       _alphaMixImage.SaveData(_currentEvent);
	if(_enableBetaMix)        _betaMixImage.SaveData(_currentEvent);
	if(_enableRBE)            _rbeImage.SaveData(_currentEvent);
	if(_enableUncertainties)  _biodoseUncertaintyImage.SaveData(_currentEvent);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void GateBioDoseActor::BeginOfRunAction(const G4Run* r) {
	GateVActor::BeginOfRunAction(r);
	GateDebugMessage("Actor", 3, "GateBioDoseActor -- Begin of Run\n");

	if(_enableEdep)     _edepImage.Reset();
	if(_enableDose)     _doseImage.Reset();
	if(_enableBioDose)  _bioDoseImage.Reset();
	if(_enableAlphaMix) _alphaMixImage.Reset();
	if(_enableBetaMix)  _betaMixImage.Reset();
	if(_enableRBE)      _rbeImage.Reset();
	if(_enableUncertainties) {
		_biodoseUncertaintyImage.Reset();
		_eventEdepImage.Reset();
		_eventAlphaImage.Reset();
		_squaredAlphaMixImage.Reset();
		_eventSqrtBetaImage.Reset();
		_squaredSqrtBetaMixImage.Reset();
		_alphaMixSqrtBetaMixImage.Reset();
	}
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void GateBioDoseActor::EndOfRunAction(const G4Run* r) {
	GateVActor::EndOfRunAction(r);
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void GateBioDoseActor::BeginOfEventAction(const G4Event* e) {
	GateVActor::BeginOfEventAction(e);
	++_currentEvent;

	_eventVoxelIndices.clear();

	_eventEdepImage.Reset();
	_eventAlphaImage.Reset();
	_eventSqrtBetaImage.Reset();
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void GateBioDoseActor::EndOfEventAction(const G4Event* e) {
	GateVActor::EndOfEventAction(e);

	for(auto const& index: _eventVoxelIndices) {
		auto const edep = _eventEdepImage.GetValue(index);
		auto const alphaMix = _eventAlphaImage.GetValue(index) / edep;
		auto const sqrtBetaMix = _eventSqrtBetaImage.GetValue(index) / edep;

		_squaredAlphaMixImage.AddValue(index, alphaMix * alphaMix);
		_squaredSqrtBetaMixImage.AddValue(index, sqrtBetaMix * sqrtBetaMix);
		_alphaMixSqrtBetaMixImage.AddValue(index, alphaMix * sqrtBetaMix);

		++_depositedMap[index].n; // TODO useless?
	}
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void GateBioDoseActor::UserSteppingActionInVoxel(const int index, const G4Step* step) {
	double const weight     = step->GetTrack()->GetWeight();
	double const energyDep  = step->GetTotalEnergyDeposit() * weight;

	if(energyDep == 0)  return;
	if(index < 0)       return;

	auto it = _depositedMap.find(index);
	if(it == std::end(_depositedMap)) {
		_depositedMap[index] = {0, 0, 0, 0, 0};
		it = _depositedMap.find(index);
	}

	auto& deposited = (*it).second;

	// Accumulate energy inconditionnaly
	deposited.energy += energyDep;

	if(_enableDose || _enableBioDose || _enableRBE) {
		decltype(_doseImage)* image = nullptr;
		if(_enableDose)         image = &_doseImage;
		else if(_enableBioDose) image = &_bioDoseImage;
		else if(_enableRBE)     image = &_rbeImage;

		auto* currentMaterial = step->GetPreStepPoint()->GetMaterial();
		double density = currentMaterial->GetDensity();
		double mass = image->GetVoxelVolume() * density;

		deposited.dose += energyDep / mass / CLHEP::gray;
	}

	// Get information from step
	// Particle
	G4int nZ = step->GetTrack()->GetDefinition()->GetAtomicNumber();
	double kineticEnergyPerNucleon = (step->GetPreStepPoint()->GetKineticEnergy()) / (step->GetTrack()->GetDefinition()->GetAtomicMass()); //OK

	++_eventCount;

	// Accumulation of alpha/beta if ion type if known
	// -> check if the ion type is known
	if(_energyMaxForZ.count(nZ) != 0) {
		++_eventWithKnownIonCount;

		//The max values in the database aren't being taking into account
		//so for now it's coded like this to be sure the code takes them into account
		double energyMax = _energyMaxForZ.at(nZ);

		AlphaBetaInterpolTable::const_iterator itr2;
		if (kineticEnergyPerNucleon >= energyMax) {
			// If the kinetic energy is the maximum value in the alpha beta tables,
			// we have to use the a and b coefficient for this maximum value
			Fragment fragmentKineticEnergyMax{nZ, energyMax};
			itr2 = _alphaBetaInterpolTable.find(fragmentKineticEnergyMax);
		} else {
			// We pair the ion type and the kinetic energy
			Fragment fragmentKineticEnergy{nZ, kineticEnergyPerNucleon};
			itr2 = _alphaBetaInterpolTable.upper_bound(fragmentKineticEnergy);
		}

		// Calculation of EZ, alphaDep and betaDep (K = a*EZ+b*E)
		auto const& interpol = (*itr2).second;

		double alpha = (interpol.alpha.a * kineticEnergyPerNucleon + interpol.alpha.b);
		double sqrtBeta = (interpol.sqrtBeta.a * kineticEnergyPerNucleon + interpol.sqrtBeta.b);

		if(alpha < 0) alpha = 0;
		if(sqrtBeta < 0) sqrtBeta = 0;

		double alphaDep = alpha * energyDep;
		double sqrtBetaDep = sqrtBeta * energyDep;

		// Accumulate alpha/beta
		deposited.alpha     += alphaDep;
		deposited.sqrtBeta  += sqrtBetaDep;

		if(_enableUncertainties) {
			_eventEdepImage.AddValue(index, energyDep);
			_eventAlphaImage.AddValue(index, alphaDep);
			_eventSqrtBetaImage.AddValue(index, sqrtBetaDep);
		}

		_eventVoxelIndices.insert(index);
	}
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void GateBioDoseActor::ResetData() {
	if(_enableEdep)     _edepImage.Reset();
	if(_enableDose)     _doseImage.Reset();
	if(_enableBioDose)  _bioDoseImage.Reset();
	if(_enableAlphaMix) _alphaMixImage.Reset();
	if(_enableBetaMix)  _betaMixImage.Reset();
	if(_enableRBE)      _rbeImage.Reset();
}
//-----------------------------------------------------------------------------
