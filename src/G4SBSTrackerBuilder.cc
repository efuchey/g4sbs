#include "G4SBSTrackerBuilder.hh"

#include "G4SBSDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4TwoVector.hh"
#include "G4GenericTrap.hh"
#include "G4Polycone.hh"
#include "G4RotationMatrix.hh"

#include "G4SBSGEMSD.hh"
#include "sbstypes.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4SBSTrackerBuilder::G4SBSTrackerBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  fFrontTrackerZoffset = 0.0;
  fFrontTrackerZoffsetIsSet = false;

  fMWPC_NplaneWires = 6;//
  fMWPC_NWiresPerPlane = 1000;//
  fMWPC_InterWireGap = 2.0*mm;//1.0*mm;//
  fMWPC_SupportThickness = 2.0*mm;//1.6*mm;//
  fMWPC_SupportWidth = 1.*cm;
  fMWPC_CathodeMylarThickness = 0.001*mm;//0.001*mm;//1.0 um Mylar
  fMWPC_CathodeAlThickness =    0.00001*mm; //0.01um Al;//
  fMWPC_WirePlaneThickness =    0.008*mm;//8 um wires;//
  fChkOvLaps = false;//true;for debugging
}

G4SBSTrackerBuilder::~G4SBSTrackerBuilder(){;}

void G4SBSTrackerBuilder::BuildComponent(G4LogicalVolume *){
  // This shouldn't get called
  ;
}

//This routine allows us to flexibly position GEM modules without code duplication:
void G4SBSTrackerBuilder::BuildComponent(G4LogicalVolume *Mother, G4RotationMatrix *rot, G4ThreeVector pos, unsigned int nplanes, vector<double> zplanes, vector<double> wplanes, vector<double> hplanes, G4String SDname, G4bool ispolarimeter, G4double yoffset ) 
{
  //This routine will create and position a GEM tracker consisting of nplanes planes centered at position pos oriented with rotation rot wrt logical volume Mother. 
  //The list of z coordinates, widths and heights of the planes are passed as arguments:

  G4bool GEMflip = fDetCon->GetFlipGEM();
  
  // How should we interpret the rotation Matrix rot? It is the rotation that orients the z axis of the tracker with respect to the mother volume. 
  // Since pos is the nominal position of the tracker with respect to the mother volume, 
  // the positioning of the centers of the tracker planes should be pos + zplane * tracker_zaxis
  G4ThreeVector zaxis(0,0,1);
  G4ThreeVector yaxis(0,1,0);
  
  zaxis *= rot->inverse();
  
  G4String TrackerPrefix = SDname; 

  //Create sensitive detector for this tracker:
  G4String GEMSDname = SDname;
  G4String GEMSDname_nopath = SDname;
  GEMSDname_nopath.erase(0,SDname.find_last_of('/')+1);
  G4String GEMcolname = GEMSDname_nopath; //We have to remove all the directory structure from the Hits collection name or else GEANT4 SDmanager routines will not handle correctly.
  GEMcolname += "HitsCollection";

  G4SBSGEMSD* GEMSD;

  if( !(GEMSD = (G4SBSGEMSD*) fDetCon->fSDman->FindSensitiveDetector(GEMSDname)) ){ //Make sure SD with this name doesn't already exist
    GEMSD = new G4SBSGEMSD( GEMSDname, GEMcolname );
    fDetCon->fSDman->AddNewDetector(GEMSD);
    (fDetCon->SDlist).insert(GEMSDname);
    fDetCon->SDtype[GEMSDname] = G4SBS::kGEM;
  }

  fDetCon->InsertSDboundaryVolume( Mother->GetName(), GEMSDname );
  
  if( !( nplanes > 0 && zplanes.size() == nplanes && wplanes.size() == nplanes && hplanes.size() == nplanes ) ){
    G4cout << "MakeTracker called with invalid arguments, doing nothing..." << G4endl;
    return;
  }

  //Best way to handle "flipping" of GEMs? Boxes are symmetric, just flip ordering of z positions!
  
  //Define z extent of various GEM layer components:

  const unsigned int nlayers = 24;

  double gempz[nlayers] = {
    120.0*um, 3.0*mm, 120.0*um, 5.0*um, 50.0*um, 3.0*mm, // cover + cathode
    5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM0
    5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM1
    5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM2
    10.0*um, 50.0*um, 180.0*um, 120.0*um, 3.0*mm, 120.0*um // readout + honeycomb
  };

  //Define GEM materials:

  G4Material *gempm[nlayers] = {
    GetMaterial("NEMAG10"), GetMaterial("NOMEX"), GetMaterial("NEMAG10"), GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("Copper"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("Copper"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("Copper"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("NEMAG10"), GetMaterial("NEMAG10"), GetMaterial("NOMEX"), GetMaterial("NEMAG10")
  };

  unsigned int gidx, gpidx;

  //Determine the largest transverse dimensions of a plane in this tracker:
  double gemmaxw = 0.0;
  double gemmaxh = 0.0;

  for( gidx = 0; gidx < nplanes; gidx++ ){
    if( wplanes[gidx] > gemmaxw ) gemmaxw = wplanes[gidx];
    if( hplanes[gidx] > gemmaxh ) gemmaxh = hplanes[gidx];
  }

  //Determine the total z extent of a plane in this tracker:
  double gempzsum = 0.0;
  for( gpidx = 0; gpidx < nlayers; gpidx++ ){
    gempzsum += gempz[gpidx];
  }

  //implement aluminum shielding: 
  G4bool useAlshield = fDetCon->GetGEMuseAlshield();

  // G4Box *alshield_box;
  // G4LogicalVolume *alshield_log;
  
  G4double Al_thick=0.0, Airgap_thick=0.0;
  if( useAlshield ){
    //G4cout << "building aluminum shielding for GEMs..."; 
    
    Al_thick = fDetCon->GetGEMAlShieldThick();
    Airgap_thick = fDetCon->GetGEMAirGapThick();
    

    //G4cout << " complete" << G4endl;
  }

  // char gpname[20][50][3][255];
  // char gemname[10][255], gemboxname[10][255], gemlogname[10][255];
  char cgidx[255];
  char cgpidx[255];

  G4Box *gpbox;
  G4LogicalVolume *gplog;

  G4VisAttributes *gemvisatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 0.0 ) );

  gemvisatt->SetForceWireframe(true);

  G4VisAttributes *gemsdvisatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 0.0 ) );

  
  
  for( gidx = 0; gidx < nplanes; gidx++ ){
    sprintf( cgidx, "%02d", gidx );

    // sprintf(gemboxname[gidx], "gembox_%02d", gidx);
    // sprintf(gemlogname[gidx], "gemlog_%02d", gidx);
    G4String gemboxname =  TrackerPrefix + G4String("_gembox_");
    gemboxname += cgidx;
    G4String gemlogname = TrackerPrefix + G4String("_gemlog_");
    gemlogname += cgidx;

    G4Box *gembox = new G4Box( gemboxname, wplanes[gidx]/2.0, hplanes[gidx]/2.0, gempzsum/2.0 );
    G4LogicalVolume *gemlog = new G4LogicalVolume( gembox, GetMaterial("Air"), gemlogname, 0, 0, 0 );

    gemlog->SetVisAttributes( gemvisatt );

    double ztemp = 0.0;
    double sign = 1.0;
    if( GEMflip ) { //this reverses ordering of GEM materials front-to-back!
      ztemp = gempzsum;
      sign = -1.0;
    }

    if( useAlshield ){ //put Aluminum shielding upstream of GEM:
      G4String shieldboxname = TrackerPrefix + G4String("alshield_gem") + cgidx + G4String("box");
      G4String shieldlogname = TrackerPrefix + G4String("alshield_gem") + cgidx + G4String("log");
      
      G4String frontshieldphysname = TrackerPrefix + G4String("alshieldfront_gem") + cgidx + G4String("phys");
      G4String backshieldphysname = TrackerPrefix + G4String("alshieldback_gem") + cgidx + G4String("phys");

      G4cout << "Building and placing GEM aluminum shielding, plane = " << gidx << G4endl;

      G4Box *alshield_box = new G4Box( shieldboxname, wplanes[gidx]/2.0, hplanes[gidx]/2.0, Al_thick/2.0 );
      G4LogicalVolume *alshield_log = new G4LogicalVolume( alshield_box, GetMaterial("Aluminum"), shieldlogname, 0, 0, 0 );

      G4VisAttributes *alshield_visatt = new G4VisAttributes( G4Colour( 0.6, 0.6, 0.6 ) );
      //alshield_visatt->SetForceWireframe( true );
      
      alshield_log->SetVisAttributes( alshield_visatt );
      
      //position the front shielding upstream of the GEM with possible air gap:
      G4ThreeVector frontshieldpos = pos + (zplanes[gidx] - gempzsum/2.0 - Airgap_thick - Al_thick/2.0) * zaxis;
      
      new G4PVPlacement( 0, frontshieldpos, alshield_log, frontshieldphysname, Mother, true, gidx+1, false );

      //position the back shielding immediately downstream of the GEM with no air gap:
      G4ThreeVector backshieldpos = frontshieldpos + (gempzsum + Al_thick + Airgap_thick ) * zaxis;

      new G4PVPlacement( 0, backshieldpos, alshield_log, backshieldphysname, Mother, true, gidx+1, false );
    }

    double zdrift; //z position of center of drift region relative to GEM plane box:
    
    for( gpidx = 0; gpidx < nlayers; gpidx++ ){
      sprintf( cgpidx, "_%02d_%03d_", gidx, gpidx );
      
      G4String gempboxname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("box");
      G4String gemplogname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("log");
      G4String gempphysname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("phy");

      ztemp += sign*gempz[gpidx]/2.0;

      gpbox = new G4Box( gempboxname, wplanes[gidx]/2.0, hplanes[gidx]/2.0, gempz[gpidx]/2.0 );
      gplog = new G4LogicalVolume( gpbox, gempm[gpidx], gemplogname, 0, 0, 0 );

      new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, ztemp - gempzsum/2.0 ), gplog, gempphysname, gemlog, false, 0, false ); 

      

      //Assign sensitive volume: why 5?  // SPR: This is the gas drift region
      if( gpidx == 5 ){
	gplog->SetSensitiveDetector(GEMSD);
	gplog->SetVisAttributes( gemsdvisatt );
	// Until we implement actual strips/wires in the GEM construction, the detmap is irrelevant for the GEMs
	zdrift = ztemp - gempzsum/2.0;
      } else {
	gplog->SetVisAttributes( G4VisAttributes::GetInvisible() );
      }

      ztemp += sign*gempz[gpidx]/2.0;
      
    }

    G4String gemname = TrackerPrefix + G4String("_gemphys_") + cgidx;
    //Now place the fully constructed GEM plane in the mother logical volume:
    //G4ThreeVector plane_pos = pos + G4ThreeVector( 0.0, 0.0, zplanes[gidx] );
    G4ThreeVector plane_pos = pos + zplanes[gidx] * zaxis + yoffset * yaxis;
    //Now we are positioning the GEM AFTER positioning all of its components inside its logical volume:
    if( gidx == 0 ){ //set Z offset of GEMSD to equal the center of the drift region of first GEM

      //NOTE: the typical wor
      
      if( ispolarimeter && fFrontTrackerZoffsetIsSet ){
	GEMSD->SetZoffset( GetFrontTrackerZoffset() );
      } else {
	GEMSD->SetZoffset( plane_pos.getZ() + zdrift );
	SetFrontTrackerZoffset( GEMSD->GetZoffset() );
      } 
    }
    new G4PVPlacement( rot, plane_pos, gemlog, gemname, Mother, true, gidx+1, false );
  }
}


//
void G4SBSTrackerBuilder::BuildComponent_MWPC(G4LogicalVolume *Mother, G4RotationMatrix *rot, G4ThreeVector pos, unsigned int nplanes, vector<double> zplanes, vector<double> wplanes, vector<double> hplanes, G4String SDname, G4bool ispolarimeter, G4double yoffset ) 
{
  //This routine will create and position a MWPC tracker consisting of nplanes planes centered at position pos oriented with rotation rot wrt logical volume Mother. 
  //The list of z coordinates, widths and heights of the planes are passed as arguments:

  // How should we interpret the rotation Matrix rot? It is the rotation that orients the z axis of the tracker with respect to the mother volume. 
  // Since pos is the nominal position of the tracker with respect to the mother volume, 
  // the positioning of the centers of the tracker planes should be pos + zplane * tracker_zaxis
  G4ThreeVector zaxis(0,0,1);
  G4ThreeVector yaxis(0,1,0);
  
  zaxis *= rot->inverse();

  double WireOffset[2] = {0.0, -fMWPC_InterWireGap*0.5};
  double WirePitch[3] = {0.0*deg, -10.0*deg, +10.0*deg};
  
  G4VisAttributes *mwpcvisatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 0.0 ) );

  mwpcvisatt->SetForceWireframe(true);

  for(int i = 0; i<nplanes; i++){
    //each wire chamber has
    G4double gap = 0.0*mm;
    G4double MWPC_z_total = (fMWPC_SupportThickness*2+fMWPC_CathodeMylarThickness+fMWPC_CathodeAlThickness*2+fMWPC_WirePlaneThickness+gap*4)*fMWPC_NplaneWires+fMWPC_CathodeMylarThickness+fMWPC_CathodeAlThickness*2+fMWPC_SupportThickness*2+gap*4;

    G4Box* MWPCmother_box = new G4Box("MWPCmother_solid",
				      wplanes[i]*0.5+fMWPC_SupportWidth+1.0*mm,
				      hplanes[i]*0.5+fMWPC_SupportWidth+1.0*mm,
				      MWPC_z_total);

    G4LogicalVolume* MWPCmother_log = 
      new G4LogicalVolume(MWPCmother_box, 
			  GetMaterial("MWPCgas"),
			  "MWPCmother_log");
    MWPCmother_log->SetVisAttributes(mwpcvisatt);
    new G4PVPlacement(0, 
		      G4ThreeVector(0.0, yoffset, zplanes[i]), 
		      MWPCmother_log,
		      "MWPCmother_phys", 
		      Mother, 
		      false, 
		      0, 
		      fChkOvLaps);
    
    G4double zpos = -MWPC_z_total*0.5 + fMWPC_SupportThickness*0.5+gap;
    G4double theta = 0.0*deg;
    G4cout << "MWPC delta Z " << MWPC_z_total << " Z position = " << zpos/mm << " mm " << endl;
    //starting with one G10 support
    BuildMWPCG10Support(MWPCmother_log, wplanes[i], hplanes[i], zpos);
    for(int k = 0; k<fMWPC_NplaneWires; k++){
      zpos+= (fMWPC_SupportThickness+fMWPC_CathodeMylarThickness)*0.5+fMWPC_CathodeAlThickness+gap;
      G4cout << " Z position = " << zpos/mm << " mm " << endl;
      BuildMWPCCathode(MWPCmother_log, wplanes[i], hplanes[i], zpos);
      zpos+= (fMWPC_CathodeMylarThickness+fMWPC_SupportThickness)*0.5+fMWPC_CathodeAlThickness+gap;
      G4cout << " Z position = " << zpos/mm << " mm " << endl;
      BuildMWPCG10Support(MWPCmother_log, wplanes[i], hplanes[i], zpos);
      zpos+= (fMWPC_SupportThickness+fMWPC_WirePlaneThickness)*0.5+gap;
      G4cout << " Z position = " << zpos/mm << " mm " << endl;
      BuildMWPCWires(MWPCmother_log, wplanes[i], hplanes[i], zpos, WirePitch[k/2], WireOffset[k%2], i*fMWPC_NplaneWires+k, SDname);
      zpos+= (fMWPC_WirePlaneThickness+fMWPC_SupportThickness)*0.5+gap;
      G4cout << " Z position = " << zpos/mm << " mm " << endl;
      BuildMWPCG10Support(MWPCmother_log, wplanes[i], hplanes[i], zpos);
    }
    //after building all wires, capping with a cathode plane and a G10 support plane
    zpos+= (fMWPC_SupportThickness+fMWPC_CathodeMylarThickness)*0.5+fMWPC_CathodeAlThickness+gap;
    G4cout << " Z position = " << zpos/mm << " mm " << endl;
    BuildMWPCCathode(MWPCmother_log, wplanes[i], hplanes[i], zpos);
    zpos+= (fMWPC_CathodeMylarThickness+fMWPC_SupportThickness)*0.5+fMWPC_CathodeAlThickness+gap;
    G4cout << " Z position = " << zpos/mm << " mm " << endl;
    BuildMWPCG10Support(MWPCmother_log, wplanes[i], hplanes[i], zpos);
  }

}



// EPAF 2026/09/16
// Methods to build MWPC design components.
void G4SBSTrackerBuilder::BuildMWPCWires(G4LogicalVolume *motherlog, G4double deltax, G4double deltay, G4double zPos, G4double angle, G4double wireoffset, G4int WirePlaneNum, G4String SDname)
{
  //fMWPC_NWiresPerPlane = deltay/fMWPC_InterWireGap-1;
  G4cout << " delta y = " << deltay/mm << " mm; number of wires per plane: " << fMWPC_NWiresPerPlane << G4endl;
  
  G4Box* MWPCWiresGas_solid = 
    new G4Box("MWPCWiresGas_solid", 
	      deltax*0.5, deltay*0.5,
	      fMWPC_SupportThickness+fMWPC_WirePlaneThickness*0.5);
  
  //mother volume for the wires...
  G4LogicalVolume* MWPCWiresGas_log = 
    new G4LogicalVolume(MWPCWiresGas_solid,
			GetMaterial("MWPCgas"),
			"MWPCWiresGas_log");
  
  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, zPos), 
		    MWPCWiresGas_log,
		    "MWPCWiresGas_phys", 
		    motherlog, 
		    false, 
		    WirePlaneNum, 
		    fChkOvLaps);

  G4VisAttributes *tpc_gas_visatt = new G4VisAttributes( G4Colour( 0.0, 0.75, 1.0) );
  tpc_gas_visatt->SetForceWireframe(true);

  G4VisAttributes *wire_visatt = new G4VisAttributes( G4Colour( 0.8, 0.75, 0.0) );
  
  MWPCWiresGas_log->SetVisAttributes( tpc_gas_visatt );
  
  //Split gas in cells for wires... will be tedious
  //TODO: introduce pitch (rn, pitch is assumed 0 degrees, which is fine for our purpose)
  G4Box* WireCell_box = 
    new G4Box("WireCell_box",
	      deltax,
	      fMWPC_InterWireGap*0.5*cos(angle),
	      fMWPC_SupportThickness+fMWPC_WirePlaneThickness*0.5);
  
  // G4Tubs* Wire_solid_base =
  //   new G4Tubs("Wire_solid_base",
  // 	       0,
  // 	       fMWPC_WirePlaneThickness*0.5,
  // 	       deltax,
  // 	       0.*deg,
  // 	       360.*deg);
  
  G4String MWPCSDname = SDname;
  G4String MWPCcolname = "MWPCHitsCollection";
  
  G4SBSGEMSD* MWPCSD;
  if( !(MWPCSD = (G4SBSGEMSD*) fDetCon->fSDman->FindSensitiveDetector(MWPCSDname)) ){ //Make sure SD with this name doesn't already exist
    MWPCSD = new G4SBSGEMSD( MWPCSDname, MWPCcolname );
    fDetCon->fSDman->AddNewDetector(MWPCSD);
    (fDetCon->SDlist).insert(MWPCSDname);
    fDetCon->SDtype[MWPCSDname] = G4SBS::kGEM;
  }
  fDetCon->InsertSDboundaryVolume( motherlog->GetName(), MWPCSDname );
 
  for(int i = 0; i<fMWPC_NWiresPerPlane; i++){
    G4double WirePosition = -fMWPC_InterWireGap*0.5*fMWPC_NWiresPerPlane+wireoffset+fMWPC_InterWireGap*i;
    //G4double WirePosition = -fMWPC_InterWireGap*0.5*1000+wireoffset+fMWPC_InterWireGap*i;
    G4cout << " wire position (mm): " << WirePosition/mm << G4endl;
    //TODO: introduce pitch (rn, pitch is assumed 0 degrees)

    G4Tubs* Wire_solid_base =
      new G4Tubs("Wire_solid_base",
		 0,
		 fMWPC_WirePlaneThickness*0.5,
		 deltax*0.5/cos(angle),
		 0.*deg,
		 360.*deg);
    
    G4RotationMatrix* RotationWire = new G4RotationMatrix();
    RotationWire->rotateY(90.0*deg);
    RotationWire->rotateX(angle);

    G4IntersectionSolid* Wire_solid = 
      new G4IntersectionSolid("Wire_solid",
			      MWPCWiresGas_solid,
			      Wire_solid_base,
			      RotationWire,
			      G4ThreeVector(0, 0, 0));
    
    G4RotationMatrix* RotationWireCell = new G4RotationMatrix();
    RotationWireCell->rotateZ(angle);
    
    G4IntersectionSolid* WireCell_solid_0 =
      new G4IntersectionSolid("WireCell_solid_0",
			      MWPCWiresGas_solid,
			      WireCell_box,
			      RotationWireCell,
			      G4ThreeVector(0, 0, 0));
    
    
    G4SubtractionSolid* WireCell_solid =
      new G4SubtractionSolid("WireCell_solid",
			     WireCell_solid_0,
			     Wire_solid_base,
			     RotationWire,
			     G4ThreeVector(0, 0, 0));
    
    G4LogicalVolume* WireCell_log = 
      new G4LogicalVolume(WireCell_solid,
			  GetMaterial("MWPCgas"),
			  "WireCell_log");
    
    G4LogicalVolume* Wire_log = 
      new G4LogicalVolume(Wire_solid_base, //Wire_solid,
			  GetMaterial("Tungsten"),
			  "Wire_log");
    
    WireCell_log->SetVisAttributes( tpc_gas_visatt );
    Wire_log->SetVisAttributes( wire_visatt );

    G4int WireUniqueID = WirePlaneNum*fMWPC_NWiresPerPlane+i;
    //   //G4cout << " j "<<  j << " WireUniqueID " << WireUniqueID << " cos angle??? " << cos(angle+j*180.0*deg) << " sin angle??? " << sin(angle+j*180.0*deg) << G4endl; 

    // G4RotationMatrix* RotationPlane = new G4RotationMatrix();
    // RotationPlane->rotateZ(-angle+j*180.0*deg);
    
    new G4PVPlacement(0, 
		      G4ThreeVector(0, WirePosition, 0),
		      WireCell_log,
		      "WireCell", 
		      MWPCWiresGas_log, 
		      false, 
		      WireUniqueID, 
		      fChkOvLaps
		      );
      
    WireCell_log->SetSensitiveDetector(MWPCSD);
    
    // RotationWire = new G4RotationMatrix();
    // RotationWire->rotateY(90.0*deg);
    // RotationWire->rotateX(angle);
      
    new G4PVPlacement(RotationWire,
		      G4ThreeVector(0, WirePosition, 0), 
		      Wire_log,
		      "Wire", 
		      MWPCWiresGas_log, 
		      false, 
		      WireUniqueID, 
		      fChkOvLaps
		      );
  }
}
//good a priori
void G4SBSTrackerBuilder::BuildMWPCCathode(G4LogicalVolume *motherlog, G4double deltax, G4double deltay, G4double zPos)
{
  G4Box* MWPCCathode_solid_base = 
    new G4Box("MWPCCathode_solid_base", 
	      deltax*0.5, deltay*0.5, 
	      fMWPC_CathodeMylarThickness*0.5+fMWPC_CathodeAlThickness);
  
  G4LogicalVolume* MWPCCathode_log = 
    new G4LogicalVolume(MWPCCathode_solid_base,
			GetMaterial("Mylar"),
			"MWPCCathode_log");

  G4Box* MWPCCathodeAlCoating_solid_base = 
    new G4Box("MWPCCathodeAlCoating_solid_base", 
	       deltax*0.5, deltay*0.5, 
	      fMWPC_CathodeAlThickness*0.5);
  
  G4LogicalVolume* MWPCCathodeAlCoating_log = 
    new G4LogicalVolume(MWPCCathodeAlCoating_solid_base,
			GetMaterial("Aluminum"),
			"MWPCCathodeAlCoating_log");
  
  for(int k = 0; k<2; k++){
    new G4PVPlacement(0, 
		      G4ThreeVector(0.0, 0.0, pow(-1, k) * (fMWPC_CathodeMylarThickness*0.5-fMWPC_CathodeAlThickness*0.5) ), 
		      MWPCCathodeAlCoating_log,
		      "MWPCCathodeAlCoating_phys", 
		      MWPCCathode_log, 
		      false, 
		      0, 
		      fChkOvLaps);
  }
  
  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, zPos), 
		    MWPCCathode_log,
		    "MWPCCathode_phys", 
		    motherlog, 
		    false, 
		    0, 
		    fChkOvLaps);

  G4VisAttributes *cathode_visatt = new G4VisAttributes( G4Colour( 0.8, 0.75, 0.0) );
  cathode_visatt->SetForceWireframe(true);
  MWPCCathode_log->SetVisAttributes( cathode_visatt );
  
  MWPCCathodeAlCoating_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
}
//good a priori
void G4SBSTrackerBuilder::BuildMWPCG10Support(G4LogicalVolume *motherlog, G4double deltax, G4double deltay, G4double zPos)
{
  G4Box* MWPCSupport_LR = 
    new G4Box("MWPCSupport_LR", 
	      fMWPC_SupportWidth*0.5, 
	      deltay*0.5+fMWPC_SupportWidth, 
	      fMWPC_SupportThickness*0.5);
  
  G4Box* MWPCSupport_TB = 
    new G4Box("MWPCSupport_TB", 
	      deltax*0.5,
	      fMWPC_SupportWidth*0.5, 
	      fMWPC_SupportThickness*0.5);
  
  
  G4LogicalVolume* MWPCSupport_LR_log = 
    new G4LogicalVolume(MWPCSupport_LR,
			GetMaterial("NEMAG10"),
			"MWPCSupport_log");
  
  G4LogicalVolume* MWPCSupport_TB_log = 
    new G4LogicalVolume(MWPCSupport_TB,
			GetMaterial("NEMAG10"),
			"MWPCSupport_log");
  
  new G4PVPlacement(0, 
		    G4ThreeVector(-(deltax+fMWPC_SupportWidth)*0.5, 0.0, zPos), 
		    MWPCSupport_LR_log,
		    "MWPCSupport_L_phys", 
		    motherlog, 
		    false, 
		    0, 
		    fChkOvLaps);

  new G4PVPlacement(0, 
		    G4ThreeVector(+(deltax+fMWPC_SupportWidth)*0.5, 0.0, zPos), 
		    MWPCSupport_LR_log,
		    "MWPCSupport_R_phys", 
		    motherlog, 
		    false, 
		    0, 
		    fChkOvLaps);

  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, -(deltay+fMWPC_SupportWidth)*0.5, zPos), 
		    MWPCSupport_TB_log,
		    "MWPCSupport_B_phys", 
		    motherlog, 
		    false, 
		    0, 
		    fChkOvLaps);

  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, +(deltay+fMWPC_SupportWidth)*0.5, zPos), 
		    MWPCSupport_TB_log,
		    "MWPCSupport_T_phys", 
		    motherlog, 
		    false, 
		    0, 
		    fChkOvLaps);

  G4VisAttributes *support_visatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 1.0) );
  support_visatt->SetForceWireframe(true);
  MWPCSupport_LR_log->SetVisAttributes( support_visatt );
  MWPCSupport_TB_log->SetVisAttributes( support_visatt );
 
}
