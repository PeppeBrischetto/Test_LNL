#ifndef TREVENTACTION_HH
#define TREVENTACTION_HH

#include <G4UserEventAction.hh>
#include <globals.hh>

class TREventAction : public G4UserEventAction
{
public:
  void EndOfEventAction(const G4Event* event) override;
  TREventAction();
  ~TREventAction(){;};
  
private:
    // Numerical IDs for hit collections (-1 means unknown yet)
  G4int fCollID_cryst {-1} ;
  G4int fCollID_cryst2 {-1};
  G4int fCollID_cryst3{-1};
  G4int fCollID_trke{-1};
  G4int fCollID_trkl{-1};
  G4int fCollID_trks{-1};
  G4int fPrintModulo;
  G4int fGoodEvents;
};

#endif
