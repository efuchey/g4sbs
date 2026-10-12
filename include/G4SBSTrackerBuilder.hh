#ifndef __G4SBSTrackerBuilder_hh
#define __G4SBSTrackerBuilder_hh


#include "G4SBSComponent.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include "sbstypes.hh"

using namespace std;

class G4SBSTrackerBuilder: public G4SBSComponent {
public:
  G4SBSTrackerBuilder(G4SBSDetectorConstruction *);
  ~G4SBSTrackerBuilder();

  void BuildComponent(G4LogicalVolume *);
  void BuildComponent(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector, unsigned int, vector<double>, vector<double>, vector<double>, G4String, G4bool ispolarimeter=false, G4double yoffset=0.0 );
  void BuildComponent_MWPC(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector, unsigned int, vector<double>, vector<double>, vector<double>, G4String, G4bool ispolarimeter=false, G4double yoffset=0.0 );
  
  inline void SetFrontTrackerZoffset( G4double z ){ fFrontTrackerZoffset = z; fFrontTrackerZoffsetIsSet = true; }
  inline G4double GetFrontTrackerZoffset() const { return fFrontTrackerZoffset; }

private:

  G4bool fChkOvLaps;
  G4double fFrontTrackerZoffset;
  G4bool fFrontTrackerZoffsetIsSet; 

  void BuildMWPCCathode(G4LogicalVolume *mother, G4double dx, G4double dy, G4double z);
  void BuildMWPCG10Support(G4LogicalVolume *mother, G4double dx, G4double dy, G4double z);
  void BuildMWPCWires(G4LogicalVolume *mother, G4double dx, G4double dy, G4double z, G4double angle, G4double wireoffset, G4int WirePlaneNum, G4String sdname);
  
  // Parameters for MWPC design
  G4int fMWPC_NplaneWires;
  G4int fMWPC_NWiresPerPlane;
  G4double fMWPC_InterWireGap;
  G4double fMWPC_SupportThickness;
  G4double fMWPC_SupportWidth;
  G4double fMWPC_CathodeMylarThickness;
  G4double fMWPC_CathodeAlThickness;
  G4double fMWPC_WirePlaneThickness;
};

#endif//__G4SBSTrackerBuilder_hh
