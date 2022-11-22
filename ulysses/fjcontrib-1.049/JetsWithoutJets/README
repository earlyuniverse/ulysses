The JetsWithoutJets package is based on the physics described in:

   Jet Observables Without Jet Algorithms
   Daniele Bertolini, Tucker Chan, and Jesse Thaler
   arXiv: 1310.7584

---

The code example_basic_usage.cc illustrates a basic usage of the package.
In particular it shows how to calculate jet multiplicity,summed scalar pt, and missing energy and their corresponding trimmed version.
It also shows how to do event trimming and jet trimming with shapes.
It can be built using

    make example_basic_usage

example_basic_usage.cc makes use of an input datafile which is provided in the
/data directory, and it should be run with

   ./example_basic_usage < ../data/single-event.dat
 
The expected output can be found in example_basic_usage.ref

A more advanced example is given in example_advanced_usage.cc, which can be built and run with 

    make example_advanced_usage
   ./example_advanced_usage < ../data/single-event.dat

and compared to the expected output example_advanced_usage.ref 
This shows how to calculate shapes with multiple pTcut and R values, axes finding, user-defined jet-like event shapes, and uses of the event-storage.
Below is a description of the main classes.


---
----------------------
Shape Trimming
----------------------

  * Shape trimming is implemented as a Selector, in analogy with other grooming procedures in FastJet.  SelectorShapeTrimming takes four arguments in its constructor.

    double Rjet;   // jet radius
    double pTcut;  // jet pT cut
    double Rsub;   // subjet radius
    double fcut;   // fractional cut on pTsub/pTjet
    Selector trimmer=SelectorShapeTrimming(Rjet,pTcut,Rsub,fcut);

  * The trimmer works by taking every particle i, and only keeping particles that satisfy the criteria pTi,Rsub > fcut * pTi,Rjet.  Here, pTi,R is the radiation contained within a radius R of particle i.
  * Alternatively, one can use trimming as a jet shape, which acts as a Transformer:

    JetShapeTrimmer jetShapeTrimmer(Rsub,fcut);

Here the pT used to apply the fcut criteria is the sum pT of the jet constituents.
----------------------
Jet-like Event Shapes
----------------------

  * The various event shapes described in the physics paper derive from JetLikeEventShape.  The included examples are:

    ShapeJetMultiplicity:  "counts" the number of jets
    ShapeScalarPt:  a.k.a. HT = sum_i pTjet_i
    ShapeScalarPtToN:  sum_i (pTjet_i)^N
    ShapeSummedMass:  Sum of invariant masses
    ShapeSummedMassSquared:  Sum of invariant mass-squared
    ShapeMissingPt:  Magnitude of missing transverse momentum
    ShapeTrimmedSubjetMultiplicity:  "counts" number of kept trimmed subjets

Other functions will be added upon request.

  * All of the JetLikeEventShapes can be called using two constructors, one that gives the default version of the event shape that only has a jet radius Rjet and pT threshold pTcut:

    ShapeJetMultiplicity Nj(Rjet,pTcut);

and one that includes trimming built-in:

    ShapeJetMultiplicity Nj_trim(Rjet,pTcut,Rsub,fcut);

Once constructed they can be used using the () or result() command acting on a vector<PseudoJet>.

  * Note that Nj_trim(inputs) is not the same as Nj(trimmer(input)).  In the case of the built-in trimming, the trimming criteria are implemented as theta functions in the event shape, so quantities like pTi,R are untrimmed.  In the case of the external trimmer, particles are removed prior to calculating pT values, so pTi,R is already trimmed.  We have left in both versions of the algorithm, since it is not clear a priori which method will be better when combined with area subtraction.  For cases like ShapeSummedMass which are quadratic in input particles, we find that applying trimming first gives better performance.
 
----------------------
Multiple pTcut values and multiple R values
----------------------

  * As described in the physics paper, it is possible to vary the pT cut or R values without having to fully recalculate the event shapes.  At the moment, only multiple pT cuts have been implemented as a general class JetLikeEventShape_MultiplePtCutValues.  After using the set_input() function on a vector<PseudoJet>, you can find the eventShapeFor() a given value of the pT cut, or the ptCutFor() a given value of the event shape.  It is also possible to get the full array of pT cut and eventShape values by calling functionArray(). functionArray() returns a vector of vector<double>, each containing a pTcut/eventShape pair.
  * For the case of jet multiplicity, there is a pre-built class for ShapeJetMultiplicity_MultiplePtCutValues.  This includes an offset value by default, such that ptCutFor(n) gives the value of the event shape at (n - 0.5), as explained in the physics paper.  This offset can be changed in the constructor.
  * The only multiple R event shape at the moment is ShapeJetMultiplicity_MultipleRValues.  For this class, there is only eventShapeFor() a given value of pT cut, because the inverse function is not single-valued.  One can still get the full array of R and eventShape values by calling functionArray().  Because of the computational costs of doing the variable R procedure, we have not implemented this as a general class.

----------------------
Axes Finding through Event Shape Density
----------------------

  * Though the event shapes themselves do not have a clustering interpretation, we can assign particles to jets using a hybrid event shape density.  The class EventShapeDensity_JetAxes implements this functionality. It takes Rjet and ptcut as inputs, and after set_inputs() is called, axes() returns the jet axes as determined by a winner-take-all recombination scheme.  Note that this axis is *not* the sum of its constituents.
  * In the future, we plan to include the option of outputting the constituents of a jet axis, such that this class behaves more similarly to a standard jet algorithm.
  * This is the least developed part of the code (and the least developed part of the theory), and is likely to evolve significantly.  Correspondingly, the user should treat these results with more caution.


----------------------
Event Storage 
----------------------

  * The EventStorage class is usually used internally to cache information about particles in the event.  Since this information is requested multiple times, caching reduces computation time. If the user wants to calculate several Jet-Like Event Shapes for a given event, the user can create an EventStorage and pass this as input to the Jet-Like Event Shapes (either to the already built-in shapes or to user-defined shapes) to reduce computation time. Event Storage should be constructed with 

    EventStorage eventStorage(Rjet,pTcut,Rsub,fcut,useLocalStorage,storeNeighbors,storeMass)
    EventStorage(Rjet,pTcut,Rsub,fcut) // using default values of bools
    EventStorage eventStorage(Rjet,pTcut,useLocalStorage,storeNeighbors,storeMass) // using default Rsub=Rjet and fcut=1.0
    EventStorage(Rjet,pTcut) // using default values of bools and Rsub=Rjet and fcut=1.0


where useLocalStorage, storeNeighbors, storeMass are bools described below.

Then the storage is established by calling

    eventStorage.establishStorage(particles)

where particles is a vector<PseudoJet>. 

  * If the user does not want to do trimming the user should use the last two constructors 
  * useLocalStorage is a bool flag that allows to choose if LocalStorage should be used.  By default it is set to true, since it almost always reduces computation time (see technical details below). 
  * storeNeighbors is a bool flag that allows to choose whether neighbors (i.e. particles within R) of each particle should be stored, by default it is set to true.
  * storeMass is a bool flag that allows to choose whether mass of the jet made of neighbors (i.e. particles within R) of each particle should be stored, by default it is set to false.
  * The function to call for passing an EventStorage as input to a JetLikeEventShape is result(eventStorage).
  * Note that a JetLikeEventShape defined with an external measurement needs neighbors to be stored in the EventStorage, while the built-in Jet-Like Event Shapes do not.

----------------------
Technical Details
----------------------

  * All of the jet-like event shapes ultimately derive from "MyFunctionOfVectorOfPseudoJets".  This class is similar to fastjet::FunctionOfPseudoJet, except it returns the result() acting on a vector<PseudoJet>.  We anticipate that this functionality will be added to FastJet in the future, hence the "My" prefix to avoid naming conflicts in the future.
  * To build one's own JetLikeEventShape, all one needs to write is a function acting on the constituents of a jet, written as a function derived from MyFunctionOfVectorOfPseudoJets<double>.  For example, FunctionUnity can be used to build ShapeJetMultiplicity.
  * To build a more complicated JetLikeEventShape, one can overload the virtual result(eventStorage).  For reasons of speed, we do not offer the user a way to avoid using the EventStorage framework.  
  * At present, JetLikeEventShape only returns a double. In the future, we may add the option for JetLikeEventShape to return an arbitrary class.
  * The WinnerTakeAllRecombiner is a recombination scheme that is valid in its own right, and could be used in ordinary jet clustering.  In the physics paper, we compare the results of EventShapeDensity_JetAxes to anti-kT clustering with WinnerTakeAllRecombiner. 
  * To improve the speed of various aspects of the program, we have a LocalStorage, which finds 2R x 2R overlapping blocks of particles to limit the number of computations that need to be done.  By default, this is always used, but can be turned off in the JetLikeEventShapes or in EventShapeDensity_JetAxes using setUseLocalStorage().  (Typically, there is no reason to turn it off, unless you mistrust our code!)  At present, it cannot be turned off for Shape Trimming, though we can add that option if users request it. 
  * We have a ParticleStorage class to cache information about individual particles. Since this information is typically requested multiple times,
caching improves speed. EventStorage is a wrapper for the set of ParticleStorages, and it manages the calculation and storing of basic quantities (such as rapidity, azimuth, pT, etc) and of quantities defined in the neighborhood of each particle (such as pt_in_Rjet, pt_in_Rsub etc..) by using the basic quantities. Also, it uses by default LocalStorage to calculate quantities defined in the neighborhoods.
  * For clarity EventStorage,ParticleStorage, and LocalStorage are defined in a separate header/source file. Those are just tools used to reduce computational time. All the core functions are defined in the main JetsWithoutJets.hh header file. 
