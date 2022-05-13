(* ::Package:: *)

(* ::Title:: *)
(*Stabilizer States Package*)


BeginPackage["StabilizerStates`"];


(* ::Section::Closed:: *)
(*General Functions*)


Swap::usage="Swap[List,FirstIndex,SecondIndex] gives back List with elements at FirstIndex and SecondIndex swapped.";
Begin["`Private`"];
Swap[OrderedList_,Element1_,Element2_]:= Module[{newlist},newlist = ReplacePart[OrderedList,{Element1-> OrderedList[[Element2]],Element2-> OrderedList[[Element1]]}];newlist];
End[];

PermutedIdentity::usage="PermutedIdentity[n,{First,Second,Third,...}] gives the permutation matrix that, when multiplying the density matrix of an n-qubit system, will reorder the qubits in the manner inputted in the second argument.";
Begin["`Private`"];
PermutedIdentity[dim_,NewOrderList_]:= Module[{OriginalList = Tuples[{0,1},dim],PermutedPiece,PermutedList={},NewIdentity = IdentityMatrix[2^dim],i,j,k},For[i=1,i< (2^dim) +1,i++,PermutedPiece = {};
For[k=1,k<Length[NewOrderList]+1,k++,AppendTo[PermutedPiece,OriginalList[[i]][[-NewOrderList[[-k]]]]]];AppendTo[PermutedList,PermutedPiece]];
If[Max[NewOrderList] <= dim,For[j=1 ,j<(2^dim) +1,j++,NewIdentity = Swap[NewIdentity,j,Position[PermutedList,OriginalList[[j]]][[1,1]]];PermutedList = Swap[PermutedList,j,Position[PermutedList,OriginalList[[j]]][[1,1]]]],Print["Invalid qubit number."],Print["Unable to interpret qubit number."]];
NewIdentity];
End[];


EndPackage[];


(* ::Section::Closed:: *)
(*Entropy Functions*)


DensityMatrix::usage="DensityMatrix[StateKet] returns the density matrix for a system in StateKet, constructed by taking the outer product of StateKet with itself.";
Begin["`Private`"];
DensityMatrix[state_]:= state.ConjugateTranspose[state];
End[];

KetProjector::usage="KetProjector[DensityMatrix] returns the ket representation of the system from DensityMatrix.";
Begin["`Private`"];
KetProjector[\[Rho]_]:=Module[{stateket={},basisbra = {ConstantArray[0,Dimensions[\[Rho]][[1]]]},indexlist = Range[Dimensions[\[Rho]][[1]]],i,j},For[i=1,i<Length[indexlist]+1,i++,AppendTo[stateket,{Sqrt[Flatten[ReplacePart[basisbra,{1,i}-> 1].\[Rho].Transpose[ReplacePart[basisbra,{1,i}-> 1]]][[1]]]}]];stateket];
End[];

TraceSystem::usage="TraceSystem[DensityMatrix,{Subsystem}] computes the reduced density matrix of Subsystem from DensityMatrix. This function was built by Mark S. Tame.";
Begin["`Private`"];
SwapParts[expr_, pos1_, pos2_] := ReplacePart[#,#, {pos1,pos2}, {pos2,pos1}]&[expr]
TraceSystem[D_, s_]:= (

Qubits=Reverse[Sort[s]];
TrkM=D;

z=(Dimensions[Qubits][[1]]+1);

For[q=1,q<z,q++,
n=Log[2,(Dimensions[TrkM][[1]])];
M=TrkM;
k=Qubits[[q]];
If[k==n,
TrkM={};
For[p=1,p<2^n+1,p=p+2,
TrkM=Append[TrkM,Take[M[[p,All]],{1,2^n,2}]+Take[M[[p+1,All]],{2,2^n,2}]];
 ],
For[j=0,j<(n-k),j++,
b={0};
For[i=1,i<2^n+1,i++,
If[(Mod[(IntegerDigits[i-1,2,n][[n]]+IntegerDigits[i-1,2,n][[n-j-1]]),2])==1 && Count[b, i]  ==0, Permut={i,(FromDigits[SwapParts[(IntegerDigits[i-1,2, n]),{n},{n-j-1}],2]+1)};
b=Append[b,(FromDigits[SwapParts[(IntegerDigits[i-1,2, n]),{n},{n-j-1}],2]+1)];
c=Range[2^n];
perm=SwapParts[c,{i},{(FromDigits[SwapParts[(IntegerDigits[i-1,2, n]),{n},{n-j-1}],2]+1)}];

M=M[[perm,perm]];

 ]    
]         ;
TrkM={};
For[p=1,p<2^n+1,p=p+2,
TrkM=Append[TrkM,Take[M[[p,All]],{1,2^n,2}]+Take[M[[p+1,All]],{2,2^n,2}]];
]
   ];
]
]
;Return[TrkM])
End[];

SvonNeumannBinary::usage="SvonNeumannBinary[StateKet,{Subsystem1,Subsystem2,...}] computes the von Neumann entropy of the requested subsystem(s) of StateKet.";
Begin["`Private`"];
SvonNeumannBinary[StateKet_,SubsystemList_]:= Module[{QubitList = Reverse[Range[Log[2,Length[StateKet]]]],TraceList,EntanglementEntropy=0,NonzeroEigenlist,i},
TraceList = Delete[QubitList,ArrayReshape[SubsystemList,{Length[SubsystemList],1}]];
NonzeroEigenlist = DeleteCases[Eigenvalues[TraceSystem[DensityMatrix[StateKet],TraceList]],0];
If[(*DiagonalizableMatrixQ[TraceSystem[DensityMatrix[StateKet],TraceList]]*)True,For[i=1,i<Length[NonzeroEigenlist]+1,i++,EntanglementEntropy = EntanglementEntropy -(NonzeroEigenlist[[i]]*Log[2,NonzeroEigenlist[[i]]])],Print["Matrix is not diagonalizable."],Print["Unknown if matrix is diagonalizable."]];EntanglementEntropy];
End[];



EntropyVectorBuilder::usage="EntropyVectorBuilder[StateKet] returns the full entropy vector for the state StateKet.";
Begin["`Private`"];
EntropyVectorBuilder[StateKet_]:=Module[{EntropyVec={},NormKet,SubsystemList = DeleteCases[Subsets[Range[Log[2,Length[StateKet]]]],{}],i},
NormKet= Transpose[{Normalize[Flatten[Transpose[StateKet]]]}];For[i=1,i< Length[SubsystemList]+1,i++,AppendTo[EntropyVec,SvonNeumannBinary[NormKet,SubsystemList[[i]]]]];ArrayReshape[EntropyVec,{Length[SubsystemList],1}]];
End[];


EntropyVectors::usage="EntropyVectors[StateSet] returns the set of full entropy vectors for all states in StateSet.";
Begin["`Private`"];
Module[{EntropyVectorList = {},k},
For[k= 1,k< Length[StateSet]+1,k++,AppendTo[EntropyVectorList,EntropyVectorBuilder[StateSet[[k]]]]];
DeleteDuplicates[EntropyVectorList]];
End[];


ReducedEntropyVectorBuilder::usage="ReducedEntropyVectorBuilder[StateKet] returns the reduced entropy vector for the state StateKet.";
Begin["`Private`"];
ReducedEntropyVectorBuilder[StateKet_]:=Module[{EntropyVec={},SubsystemList = Drop[DeleteCases[Subsets[Range[Log[2,Length[StateKet]]]],{}],-((2^(Log[2,Length[StateKet]])-1)+1)/2],i},

For[i=1,i< Length[SubsystemList]+1,i++,AppendTo[EntropyVec,{SvonNeumannBinary[Transpose[{Normalize[Flatten[Transpose[StateKet]]]}],SubsystemList[[i]]]}]];
EntropyVec];
End[];


ReducedEntropyVectors::usage="ReducedEntropyVectors[StateSet] returns the set of reduced entropy vectors for all states in StateSet.";
Begin["`Private`"];
ReducedEntropyVectors[StateSet_]:= Module[{ReducedEntropyVectorList = {},k},
For[k= 1,k< Length[StateSet]+1,k++,AppendTo[ReducedEntropyVectorList,ReducedEntropyVectorBuilder[StateSet[[k]]]]];
ReducedEntropyVectorList];
End[];


StateCounterByEntropy::usage="StateCounterByEntropy[StateSet,DesiredEntropyVector] returns a count of states in StateSet with the entropy vector DesiredEntropyVector.";
Begin["`Private`"];
StateCounterByEntropy[SetOfStates_,DesiredEntropyVector_]:=Module[{CountList = {},i},
For[i=1,i<Length[SetOfStates]+1,i++,If[EntropyVectorBuilder[SetOfStates[[i]]]== DesiredEntropyVector,AppendTo[CountList,SetOfStates[[i]]],Null,Print["Unable to interpret entropy vector."]]];
CountList];
End[];



(* ::Section::Closed:: *)
(*Quantum Logic Gates*)


HadamardGate::usage= "HadamardGate[StateKet, n] returns a ket for the system StateKet after a Hadamard gate is applied to the nth qubit of StateKet.";
Begin["`Private`"];
HadamardGate[StateKet_,AffectedQubit_]:=Module[{},HadamardAction[Log[2,Length[StateKet]],AffectedQubit].StateKet];
End[];

NormalizedHadamardGate::usage= "NormalizedHadamardGate[StateKet, n] returns a normalized ket for the system StateKet after a Hadamard gate is applied to the nth qubit of StateKet.";
Begin["`Private`"];
NormalizedHadamardGate[StateKet_,AffectedQubit_]:=Module[{},Transpose[{Normalize[Flatten[Transpose[HadamardAction[Log[2,Length[StateKet]],AffectedQubit].StateKet]]]}]];
End[];

HadamardAction::usage= "HadamardAction[dim, n] returns the matrix that produces Hadamard action on the nth qubit in a dim-qubit system.";
Begin["`Private`"];
HadamardAction[Dimension_,AffectedQubit_]:=Module[{PermIdentity,BaseHMatrix = KroneckerProduct[IdentityMatrix[2^(Dimension-1)],Sqrt[2]*HadamardMatrix[2]],OriginalOrder = Range[Dimension],i},
PermIdentity = PermutedIdentity[Dimension,Swap[OriginalOrder,1,AffectedQubit]];
ConjugateTranspose[PermIdentity].BaseHMatrix.PermIdentity];
End[];

UnequalHadamardGate::usage= "UnequalHadamardGate[StateKet,n,P] returns a ket for the system StateKet after a Hadamard gate with first diagonal element P, and second diagonal element (1-P) is applied to the nth qubit of StateKet.";
Begin["`Private`"];
UnequalHadamardGate[StateKet_,AffectedQubit_,Probability1_]:=Module[{FinalState,PermIdentity,BaseHMatrix = {{Sqrt[Probability1],Sqrt[1-Probability1]},{Sqrt[1-Probability1],-Sqrt[Probability1]}},OriginalOrder = Range[Log[2,Length[StateKet]]],i},
For[i=1,i< Log[2,Length[StateKet]],i++,BaseHMatrix = KroneckerProduct[IdentityMatrix[2],BaseHMatrix]];
PermIdentity = PermutedIdentity[Log[2,Length[StateKet]],Swap[OriginalOrder,1,AffectedQubit]];
FinalState = ConjugateTranspose[PermIdentity].(BaseHMatrix.(PermIdentity.StateKet));
FinalState];
End[];

CNOTGate::usage= "CNOTGate[StateKet,c,t] returns a ket for the system StateKet after a CNOT gate with qubit number c the control bit and qubit number t the target bit.";
Begin["`Private`"];
CNOTGate[StateKet_,controlbitnumber_,targetbitnumber_]:= Module[{FinalState,PermIdentity,BaseCNOTMatrix = {{1,0,0,0},{0,0,0,1},{0,0,1,0},{0,1,0,0}},OriginalOrder = Range[Log[2,Length[StateKet]]],i},
For[i=1,i< Log[2,Length[StateKet]]-1,i++,BaseCNOTMatrix = KroneckerProduct[IdentityMatrix[2],BaseCNOTMatrix]];
PermIdentity = PermutedIdentity[Log[2,Length[StateKet]],Join[{controlbitnumber,targetbitnumber},Delete[OriginalOrder,{{controlbitnumber},{targetbitnumber}}]]];
FinalState = ConjugateTranspose[PermIdentity].(BaseCNOTMatrix.(PermIdentity.StateKet));
FinalState];
End[];

CNOTAction::usage= "CNOTAction[dim, {c,t}] returns the matrix that produces CNOT action on the qubit pair {c,t} in a dim-qubit system.";
Begin["`Private`"];
CNOTAction[Dimension_,CTpair_]:=Module[{FinalState,PermIdentity,BaseCNOTMatrix = KroneckerProduct[IdentityMatrix[2^(Dimension-2)],{{1,0,0,0},{0,0,0,1},{0,0,1,0},{0,1,0,0}}],OriginalOrder = Range[Dimension],i},
PermIdentity = PermutedIdentity[Dimension,Join[CTpair,Delete[OriginalOrder,{{CTpair[[1]]},{CTpair[[2]]}}]]];
ConjugateTranspose[PermIdentity].BaseCNOTMatrix.PermIdentity];
End[];

PhaseGate::usage= "PhaseGate[StateKet, n] returns a ket for the system StateKet after a phase gate is applied to the nth qubit of StateKet.";
Begin["`Private`"];
PhaseGate[StateKet_,AffectedQubit_]:=Module[{},PhaseAction[Log[2,Length[StateKet]],AffectedQubit].StateKet];
End[];

PhaseAction::usage= "PhaseAction[dim, n] returns the matrix that produces phase action on the nth qubit in a dim-qubit system.";
Begin["`Private`"];
PhaseAction[Dimension_,AffectedQubit_]:=Module[{PermIdentity,BasePMatrix = KroneckerProduct[IdentityMatrix[2^(Dimension-1)],{{1,0},{0,I}}],OriginalOrder = Range[Dimension],i},
PermIdentity = PermutedIdentity[Dimension,Swap[OriginalOrder,1,AffectedQubit]];
ConjugateTranspose[PermIdentity].BasePMatrix.PermIdentity];
End[];

ControlledPhaseGate::usage= "ControlledPhaseGate[StateKet, n] returns a ket for the system StateKet after a controlled-phase gate is applied to the nth qubit of StateKet.";
Begin["`Private`"];
ControlledPhaseGate[StateKet_,controlbitnumber_,targetbitnumber_]:= Module[{FinalState,PermIdentity,BaseCNOTMatrix = {{1,0,0,0},{0,0,1,0},{0,0,0,-1},{0,1,0,0}},OriginalOrder = Range[Log[2,Length[StateKet]]],i},
For[i=1,i< Log[2,Length[StateKet]]-1,i++,BaseCNOTMatrix = KroneckerProduct[IdentityMatrix[2],BaseCNOTMatrix]];
PermIdentity = PermutedIdentity[Log[2,Length[StateKet]],Join[{controlbitnumber,targetbitnumber},Delete[OriginalOrder,{{controlbitnumber},{targetbitnumber}}]]];
FinalState = ConjugateTranspose[PermIdentity].(BaseCNOTMatrix.(PermIdentity.StateKet));
FinalState];
End[];


(* ::Section::Closed:: *)
(*Quantum Circuits (All stab. circuit functions into 1, figure out what GeneratorSetVacuum does.)*)


StabilizerCircuit::usage= "StabilizerCircuit[StateKet, {OperationsList}, {QubitList}] evolves the system StateKet through a quantum ciruit of gates {OperationsList} acting on respective qubits {QubitList}. {OperationsList} is time-ordered list of gates to be applied with the first list element corresponding to the first action taken, e.g. {HadamardGate,CNOTGate,HadamardGate,PhaseGate,...}. {QubitList} is an order-list of qubit indices, including pairs, on which the gates acts in {OperationsList} act, e.g. {2,{1,2},3,2,...}.";
Begin["`Private`"];
StabilizerCircuit[StateKet_,OperationsList_,QubitList_]:= Module[{UpdatedKet = StateKet,i,t},For[t=1,t<Length[OperationsList]+1,t++,If[OperationsList[[t]]== CNOTGate,UpdatedKet = OperationsList[[t]][UpdatedKet,QubitList[[t,1]],QubitList[[t,2]]],UpdatedKet =OperationsList[[t]][UpdatedKet,QubitList[[t,1]]],UpdatedKet =OperationsList[[t]][UpdatedKet,QubitList[[t,1]]]];
Print["At time t= ",t," the system is in the state ",UpdatedKet//MatrixForm, " with entropy vector ",EntropyVectorBuilder[UpdatedKet]//MatrixForm]];
UpdatedKet];
End[];

StabilizerCircuitOnlyVector::usage= "StabilizerCircuitOnlyVector[StateKet, OperationsList,QubitList] returns a ket vector representing the system after a Hadamard gate is applied to the nth qubit of StateKet.";
Begin["`Private`"];
StabilizerCircuitOnlyVector[StateKet_,OperationsList_,QubitList_]:= Module[{UpdatedKet = StateKet,i,t},For[t=1,t<Length[OperationsList]+1,t++,If[OperationsList[[t]]== CNOTGate,UpdatedKet = OperationsList[[t]][UpdatedKet,QubitList[[t,1]],QubitList[[t,2]]],UpdatedKet =OperationsList[[t]][UpdatedKet,QubitList[[t,1]]],UpdatedKet =OperationsList[[t]][UpdatedKet,QubitList[[t,1]]]]];
((UpdatedKet[[FirstPosition[Flatten[Transpose[UpdatedKet]],_?(##!= 0&)][[1]]]][[1]])^(-1))*UpdatedKet];
End[];

StabilizerCircuitIndices::usage= "StabilizerCircuitIndices[StateKet, n] returns a ket vector representing the system after a Hadamard gate is applied to the nth qubit of StateKet.";
Begin["`Private`"];
StabilizerCircuitIndices[StateKet_,OperationsList_,QubitList_,StabilizerStateSet_]:= Module[{UpdatedKet = StateKet,IndexSet = {Position[StabilizerStateSet,StateKet]},i,t},
For[t=1,t<Length[OperationsList]+1,t++,
If[OperationsList[[t]]== CNOTGate,
UpdatedKet = OperationsList[[t]][UpdatedKet,QubitList[[t,1]],QubitList[[t,2]]],UpdatedKet =OperationsList[[t]][UpdatedKet,QubitList[[t,1]]],UpdatedKet =OperationsList[[t]][UpdatedKet,QubitList[[t,1]]]];
AppendTo[IndexSet,Position[StabilizerStateSet,((UpdatedKet[[FirstPosition[Flatten[Transpose[UpdatedKet]],_?(##!= 0&)][[1]]]][[1]])^(-1))*UpdatedKet]]];
Flatten[IndexSet]];
End[];

GeneratorSetVacuum::usage= "GeneratorSetVacuum[StateKet, n] returns a ket vector representing the system after a Hadamard gate is applied to the nth qubit of StateKet.";
Begin["`Private`"];
GeneratorSetVacuum[Dimensions_,OperationsList_,QubitList_ ]:=Module[{InitialTable={},StabilizerSet={},i ,j,k},For[i=1,i<(Dimensions+1),i++,AppendTo[InitialTable,ReplacePart[ConstantArray[0,2*Dimensions],(Dimensions+i)-> 1]]];
For[j=1,j< (Length[OperationsList]+1),j++,If[TrueQ[OperationsList[[j]]== HadamardGate],InitialTable[[All,{QubitList[[j,1]],Dimensions+QubitList[[j,1]]}]] = InitialTable[[All,{Dimensions+QubitList[[j,1]],QubitList[[j,1]]}]],
If[TrueQ[OperationsList[[j]]== CNOTGate],InitialTable[[All,QubitList[[j,2]]]]=BitXor[InitialTable[[All,QubitList[[j,1]]]],InitialTable[[All,QubitList[[j,2]]]]];InitialTable[[All,Dimensions + QubitList[[j,1]]]]=BitXor[InitialTable[[All,Dimensions+QubitList[[j,2]]]],InitialTable[[All,Dimensions + QubitList[[j,1]]]]],
If[TrueQ[OperationsList[[j]]== PhaseGate],InitialTable[[All,Dimensions + QubitList[[j,1]]]]=BitXor[InitialTable[[All,QubitList[[j,1]]]],InitialTable[[All,Dimensions+QubitList[[j,1]]]]]]]]];
Print[InitialTable//MatrixForm];
For[k=1,k<(Dimensions+1),k++,AppendTo[StabilizerSet,Replace[Replace[InitialTable[[k]],1->PauliMatrix[3],{1} ],0-> IdentityMatrix[2],{1}][[Dimensions+1;;2*Dimensions]]]];
Print["Generating set is ", StabilizerSet]];
End[];


CircuitBuilder::usage= "CircuitBuilder[StateKet, n] returns a ket vector representing the system after a Hadamard gate is applied to the nth qubit of StateKet.";
Begin["`Private`"];
CircuitBuilder[StartingState_,ChosenSubgraph_,ChosenCircuit_]:=Module[{CircuitList = {StartingState},CircuitIndex,TestState},For[CircuitIndex = 1,CircuitIndex<Length[ChosenCircuit]+1,CircuitIndex++,TestState = Apply[Dot,Reverse[ThreeQubitHamiltonianCircuit[[1;;CircuitIndex]]]].StartingState;AppendTo[CircuitList,((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState]];
CircuitList];
End[];




(* ::Section:: *)
(*Stabilizer State Functions (Combine some elements of set generators)*)


StabilizerStateSetGenerator::usage= "StabilizerStateSetGenerator[dim] returns the full dim-qubit set of stabilizer states. Set is returned as a list of vectors normalized to the first component.";
Begin["`Private`"];
StabilizerStateSetGenerator[Dimension_]:=Block[{PredictedCount = ((2^Dimension)*Product[(2^num)+1,{num,1,Dimension}]),UpdatedStabilizerStateSet={Transpose[{ReplacePart[ConstantArray[0,2^Dimension],1->1]}]},CheckSet={Transpose[{ReplacePart[ConstantArray[0,2^Dimension],1->1]}]},NormSet={},TestState ,SingleQubits = Permutations[Range[Dimension],{1}],PairedQubits = Permutations[Range[Dimension],{2}],TotalGateSet={},OperationSet={},ii,jj,kk,ll},

For[ii=1,ii<Length[SingleQubits]+1,ii++,AppendTo[TotalGateSet,HadamardAction[Dimension,SingleQubits[[ii,1]]]];AppendTo[TotalGateSet,PhaseAction[Dimension,SingleQubits[[ii,1]]]]];

For[jj=1,jj<Length[PairedQubits]+1,jj++,AppendTo[TotalGateSet,CNOTAction[Dimension,PairedQubits[[jj]]]]];

While[Length[UpdatedStabilizerStateSet]<PredictedCount,

For[kk=1,kk<Length[TotalGateSet]+1,kk++,
For[ll=1,ll<Length[CheckSet]+1,ll++,TestState = TotalGateSet[[kk]].CheckSet[[ll]];AppendTo[NormSet,((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState]]];

CheckSet = DeleteDuplicates[NormSet];

UpdatedStabilizerStateSet = DeleteDuplicates[Join[UpdatedStabilizerStateSet,CheckSet]];
Print[Length[UpdatedStabilizerStateSet]];
NormSet={};
];
UpdatedStabilizerStateSet];
End[];




NumberSelectedStateGenerator::usage= "NumberSelectedStateGenerator[StartingState, NumberDesired] returns a list of the first NumberDesired states from StartingState.";
Begin["`Private`"];
NumberSelectedStateGenerator[StartingAddress_,NumberOfDesiredStates_]:=Block[{UpdatedStabilizerStateSet={StartingAddress},CheckSet={StartingAddress},NormSet={},TestState ,SingleQubits = Permutations[Range[Log[2,Length[StartingAddress]]],{1}],PairedQubits = Permutations[Range[Log[2,Length[StartingAddress]]],{2}],TotalGateSet={},ii,jj,kk,ll},

For[ii=1,ii<Length[SingleQubits]+1,ii++,AppendTo[TotalGateSet,HadamardAction[Log[2,Length[StartingAddress]],SingleQubits[[ii,1]]]];AppendTo[TotalGateSet,PhaseAction[Log[2,Length[StartingAddress]],SingleQubits[[ii,1]]]]];

For[jj=1,jj<Length[PairedQubits]+1,jj++,AppendTo[TotalGateSet,CNOTAction[Log[2,Length[StartingAddress]],PairedQubits[[jj]]]]];

While[Length[UpdatedStabilizerStateSet]<NumberOfDesiredStates,

For[kk=1,kk<Length[TotalGateSet]+1,kk++,
For[ll=1,ll<Length[CheckSet]+1,ll++,TestState = TotalGateSet[[kk]].CheckSet[[ll]];AppendTo[NormSet,((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState]]];

CheckSet = DeleteDuplicates[NormSet];

UpdatedStabilizerStateSet = DeleteDuplicates[Join[UpdatedStabilizerStateSet,CheckSet]];
Print[Length[UpdatedStabilizerStateSet]];
NormSet={};
];
UpdatedStabilizerStateSet];
End[];


NumberSelectedStateGeneratorWrite::usage= "NumberSelectedStateGeneratorWrite[StartingAddress_,NumberOfDesiredStates_,FileDestination_] documentation coming soon...";
Begin["`Private`"];
NumberSelectedStateGeneratorWrite[StartingAddress_,NumberOfDesiredStates_,FileDestination_]:=Block[{UpdatedStabilizerStateSet={StartingAddress},CheckSet={StartingAddress},NormSet={},TestState ,SingleQubits = Permutations[Range[Log[2,Length[StartingAddress]]],{1}],PairedQubits = Permutations[Range[Log[2,Length[StartingAddress]]],{2}],TotalGateSet={},ii,jj,kk,ll},

For[ii=1,ii<Length[SingleQubits]+1,ii++,AppendTo[TotalGateSet,HadamardAction[Log[2,Length[StartingAddress]],SingleQubits[[ii,1]]]];AppendTo[TotalGateSet,PhaseAction[Log[2,Length[StartingAddress]],SingleQubits[[ii,1]]]]];

For[jj=1,jj<Length[PairedQubits]+1,jj++,AppendTo[TotalGateSet,CNOTAction[Log[2,Length[StartingAddress]],PairedQubits[[jj]]]]];

While[Length[UpdatedStabilizerStateSet]<NumberOfDesiredStates,

For[kk=1,kk<Length[TotalGateSet]+1,kk++,
For[ll=1,ll<Length[CheckSet]+1,ll++,TestState = TotalGateSet[[kk]].CheckSet[[ll]];AppendTo[NormSet,((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState]]];

CheckSet = DeleteDuplicates[NormSet];

UpdatedStabilizerStateSet = DeleteDuplicates[Join[UpdatedStabilizerStateSet,CheckSet]];
Print[Length[UpdatedStabilizerStateSet]];
NormSet={};
];
Put[UpdatedStabilizerStateSet,FileDestination];
UpdatedStabilizerStateSet];
End[];




GateSelectedStateGenerator::usage= "GGateSelectedStateGenerator[StartingAddress_,DesiredGateSet_,NumberOfDesiredStates_] documentation coming soon...";
Begin["`Private`"];
GateSelectedStateGenerator[StartingAddress_,DesiredGateSet_,NumberOfDesiredStates_]:=Block[{UpdatedStabilizerStateSet={StartingAddress},CheckSet={StartingAddress},NormSet={},TestState ,SingleQubits = Permutations[Range[Log[2,Length[StartingAddress]]],{1}],PairedQubits = Permutations[Range[Log[2,Length[StartingAddress]]],{2}],TotalGateSet={},BooleanCheck=False,ii,jj,kk,ll},

For[kk=1,kk<Length[DesiredGateSet]+1,kk++,If[DesiredGateSet[[kk,1]]=="H",AppendTo[TotalGateSet,HadamardAction[(Log[2,Length[StartingAddress]]),DesiredGateSet[[kk,2,1]]]],If[DesiredGateSet[[kk,1]]=="P",AppendTo[TotalGateSet,PhaseAction[(Log[2,Length[StartingAddress]]),DesiredGateSet[[kk,2,1]]]],AppendTo[TotalGateSet,CNOTAction[(Log[2,Length[StartingAddress]]),DesiredGateSet[[kk,2]]]],Print["Unable to interpret gate."]]]];

While[And[Length[UpdatedStabilizerStateSet]<NumberOfDesiredStates,BooleanCheck== False],

For[kk=1,kk<Length[TotalGateSet]+1,kk++,
For[ll=1,ll<Length[CheckSet]+1,ll++,TestState = TotalGateSet[[kk]].CheckSet[[ll]];AppendTo[NormSet,((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState]]];

CheckSet = DeleteDuplicates[NormSet];

UpdatedStabilizerStateSet = DeleteDuplicates[Join[UpdatedStabilizerStateSet,CheckSet]];
NormSet={};
BooleanCheck = SubsetQ[UpdatedStabilizerStateSet,CheckSet];
];
UpdatedStabilizerStateSet];
End[];

StabilizerStateGeneratorWithEntropy::usage= "StabilizerStateGeneratorWithEntropy[SetOfStates_]";
Begin["`Private`"];
StabilizerStateGeneratorWithEntropy[SetOfStates_]:=Block[{dim = Log[2,Length[SetOfStates[[1]]]],SingleQubits = Permutations[Range[Log[2,Length[SetOfStates[[1]]]]],{1}],PairedQubits = Permutations[Range[Log[2,Length[SetOfStates[[1]]]]],{2}],TotalGateSet={},TestState,StatesWithGates = {},ii,jj,kk,ll,mm,LoopCounter = 0},

For[ii=1,ii<Length[SingleQubits]+1,ii++,AppendTo[TotalGateSet,HadamardAction[dim,SingleQubits[[ii,1]]]]];

For[jj=1,jj<Length[SingleQubits]+1,jj++,AppendTo[TotalGateSet,PhaseAction[dim,SingleQubits[[jj,1]]]]];

For[mm=1,mm<Length[PairedQubits]+1,mm++,AppendTo[TotalGateSet,CNOTAction[dim,PairedQubits[[mm]]]]];

Print["Gates Compiled."];

For[kk=1,kk< Length[SetOfStates]+1,kk++,

ll=1;

While[ll<Log[2,Length[SetOfStates[[1]]]]+1,TestState = TotalGateSet[[ll]].SetOfStates[[kk]];
AppendTo[StatesWithGates,{((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState,{"H",{ll}},SetOfStates[[kk]],ReducedEntropyVectorBuilder[TestState]}];ll++];

While[ll<(2*(Log[2,Length[SetOfStates[[1]]]]))+1,TestState = TotalGateSet[[ll]].SetOfStates[[kk]];
AppendTo[StatesWithGates,{((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState,{"P",{ll-dim}},SetOfStates[[kk]],ReducedEntropyVectorBuilder[TestState]}];ll++];

While[ll<Length[TotalGateSet]+1,TestState = TotalGateSet[[ll]].SetOfStates[[kk]];
AppendTo[StatesWithGates,{((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState,{"CNOT",PairedQubits[[ll-(2*dim)]]},SetOfStates[[kk]],ReducedEntropyVectorBuilder[TestState]}];ll++];

If[LoopCounter>1000,Print[kk];LoopCounter = 0,LoopCounter +=1,Print["Confused."]];
];
StatesWithGates];
End[];


GroupBuilder::usage="GroupBuilder[List,FirstIndex,SecondIndex]";
Begin["`Private`"];
GroupBuilder[OperationList_]:=Module[{i,j,TestGroup = {},OldSet={},FinalGroup= OperationList},

While[SubsetQ[TestGroup,FinalGroup]== False,

TestGroup = DeleteDuplicates[Join[TestGroup,FinalGroup]];

For[i=1,i<Length[TestGroup]+1,i++,
For[j=1,j<Length[TestGroup]+1,j++,
AppendTo[OldSet,TestGroup[[i]].TestGroup[[j]]]]];

FinalGroup=DeleteDuplicates[OldSet];
OldSet={};
];
DeleteDuplicates[TestGroup ]];
End[];


StabilizingSubgroupFinder::usage="StabilizingSubgroupFinder[Group_,State_]";
Begin["`Private`"];
StabilizingSubgroupFinder[Group_,State_]:=Module[{StabilizerGroup={},TestState,i},
For[i=1,i<Length[Group]+1,i++,
TestState=(Group[[i]]).State;
If[(((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState) == State,AppendTo[StabilizerGroup,Group[[i]]],Null]];
StabilizerGroup];
End[];


LoopCountWords::usage="LoopCountWords[StatesWithGates_,GateSet_]";
Begin["`Private`"];
LoopCountWords[StatesWithGates_,GateSet_]:=Module[{Loops={},i,j},
For[i=1,i<Length[StatesWithGates]+1,i++, If[StatesWithGates[[i,1]]== StatesWithGates[[i,3]],AppendTo[Loops,StatesWithGates[[i]]],Null]];
Print["There are ", Length[Loops], " total loops."];

For[j=1,j<Length[GateSet]+1,j++,Print[Count[Loops,GateSet[[j]],2]," ", ToString[GateSet[[j]]], " loops."]];
];
End[];


LoopCountLoops::usage="LoopCountLoops[StatesWithGates_,GateSet_]";
Begin["`Private`"];
LoopCountLoops[StatesWithGates_,GateSet_]:=Module[{Loops={},i,j},
For[i=1,i<Length[StatesWithGates]+1,i++, If[StatesWithGates[[i,1]]== StatesWithGates[[i,3]],AppendTo[Loops,StatesWithGates[[i]]],Null]];
Loops
];
End[];



StabilizingSubgroupFinder::usage="StabilizingSubgroupFinder[Group_,State_]";
Begin["`Private`"];
StabilizingSubgroupFinder[Group_,State_]:=Module[{StabilizerGroup={},TestState,i},
For[i=1,i<Length[Group]+1,i++,
TestState=(Group[[i]]).State;
If[(((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState) == State,AppendTo[StabilizerGroup,Group[[i]]],Null]];
StabilizerGroup];
End[];


(* ::Section:: *)
(*Visualization Functions*)


(* ::Subsection:: *)
(*Graph Builders*)


SimpleReachabilityGraphBuilder::usage= "SimpleReachabilityGraphBuilder[SetofStates_,SetOfStatesWithGates_,GatesToBeDisplayed_]";
Begin["`Private`"];
SimpleReachabilityGraphBuilder[SetofStates_,SetOfStatesWithGates_,GatesToBeDisplayed_]:=Module[{TestValue,OperationsConsidered={},VertexSet = Range[Length[SetofStates]],EdgeSet={},i,j,k=0,m,G},

For[i=1,i<Length[GatesToBeDisplayed]+1,i++,
AppendTo[OperationsConsidered,Position[SetOfStatesWithGates,GatesToBeDisplayed[[i]]]]];

For[j=1,j<Length[Flatten[OperationsConsidered,1]]+1,j++,
TestValue=SetOfStatesWithGates[[Flatten[OperationsConsidered,1][[j]][[1]]]];

If[TestValue[[2,1]]=="P" ,AppendTo[EdgeSet,DirectedEdge[Position[SetofStates,TestValue[[3]]][[1]],Position[SetofStates,TestValue[[1]]][[1]]]],
AppendTo[EdgeSet,Sort[Position[SetofStates,TestValue[[3]]][[1]]<-> Position[SetofStates,TestValue[[1]]][[1]]]],Print["Cannot determine gate input."]];

k++;

If[k>1000,Print[j," of ",Length[Flatten[OperationsConsidered,1]]," edges added."];k=0,Null]];

(*Uses the edge set compiled above to render the graph.*)
G = SimpleGraph[EdgeSet];
G];
End[];

ReachabilityGraphBuilder::usage= "ReachabilityGraphBuilder[SetOfStatesWithGates_,GatesToBeDisplayed_,EntropyVectorList_,EdgeLabelBoolean_,VertexLabelChoice_,VertexCoords_]";
Begin["`Private`"];
ReachabilityGraphBuilder[SetOfStatesWithGates_,GatesToBeDisplayed_,EntropyVectorList_,EdgeLabelBoolean_,VertexLabelChoice_,VertexCoords_]:=Block[{TestValue,dimension = Length[SetOfStatesWithGates[[1,1]]],BasisKets=Tuples[{0,1},Length[SetOfStatesWithGates[[1,1]]]],OperationsConsidered={},EdgeList={},EdgeStyleList = {Dashing[{}],Dashing[{Small}],Dashing[{Tiny}]},EdgeLabelSet={},VList={},VertexStyleSet={},VertexColorChoices = {Blue,Red,Orange,Purple,Yellow,Black,Cyan,Pink,Brown,LightGreen,White,LightRed,LightBlue,LightYellow,LightMagenta,LightOrange,LightGray,LightCyan,LightBrown},VertexLabelSet = None,h,i,j,k,l,m,LoopCounter=0,G},

For[h=1,h<Length[BasisKets]+1,h++,ReplacePart[BasisKets,h-> Ket[BasisKets[[h]]]]];

(*All the relevant operations, those represented by edges in the graph, are compiled.*)
For[i=1,i<Length[GatesToBeDisplayed]+1,i++,
OperationsConsidered = Join[OperationsConsidered,Most/@Position[SetOfStatesWithGates,GatesToBeDisplayed[[i]]]]];

(*Operation list is reformatted and edges generated/assigned.*)
For[j=1,j<Length[OperationsConsidered]+1,j++,
TestValue=SetOfStatesWithGates[[OperationsConsidered[[j,1]]]];

AppendTo[VList,TestValue[[1]]];

If[TestValue[[2,1]]=="P" ,AppendTo[EdgeList,Style[DirectedEdge[TestValue[[3]],TestValue[[1]]],{EdgeStyleList[[TestValue[[2,2]]]]}]],
If[TestValue[[2,1]]=="H",AppendTo[EdgeList,Style[Sort[TestValue[[3]]<-> TestValue[[1]]],{Green,EdgeStyleList[[TestValue[[2,2]]]]}]],AppendTo[EdgeList,Style[Sort[TestValue[[3]]<-> TestValue[[1]]],{Magenta,EdgeStyleList[[TestValue[[2,2]]]]}]]]];

AppendTo[EdgeLabelSet,EdgeList[[j]]-> ToString[TestValue[[2]]]];

AppendTo[VertexStyleSet,TestValue[[1]] -> VertexColorChoices[[Position[EntropyVectorList,TestValue[[4]]][[1,1]]]]];

LoopCounter++;

(*For user sanity.*)
If[LoopCounter>1000,Print[j," of ",Length[OperationsConsidered]," edges added."];LoopCounter=0,Null]];

(*Generate vertex labels desired.*)
If[VertexLabelChoice== "Ket",
VList=DeleteDuplicates[VList];
VertexLabelSet = {};
For[k=1,k<Length[VList]+1,k++,
For[m=1,m<dimension+1,m++,AppendTo[VertexLabelSet,ReplacePart[VList[[k]],m-> Ket[VList[[k,m]]]]]]],
If[VertexLabelChoice== "Address",
VList=DeleteDuplicates[VList];
VertexLabelSet = {};
For[k=1,k<Length[VList]+1,k++,
AppendTo[VertexLabelSet,VList[[k]]-> ToString[Flatten[Transpose[VList[[k]]]]]]],Null]];

(*Generate edge labels if desired.*)
If[EdgeLabelBoolean== True,
EdgeLabelSet = DeleteDuplicates[EdgeLabelSet],
EdgeLabelSet = None];

Graph[Sort[DeleteDuplicates[VList]],DeleteDuplicates[EdgeList],VertexStyle-> Sort[DeleteDuplicates[VertexStyleSet]],VertexLabels->VertexLabelSet,VertexLabelStyle->Directive[Italic,Thick,7],VertexCoordinates->VertexCoords,EdgeLabels->EdgeLabelSet]
];
End[];

OperationLocator::usage= "OperationLocator[StatesAndGates_,OperationList_]";
Begin["`Private`"];
OperationLocator[StatesAndGates_,OperationList_]:=Module[{PositionList={},i,j},
For[i=1,i<Length[OperationList]+1,i++,PositionList = Join[PositionList,Most/@Position[StatesAndGates,OperationList[[i]]]]];
StatesAndGates[[Sort[Flatten[PositionList]]]]];
End[];

GraphFinderExactVertex::usage= "GraphFinderExactVertex[InitialState_,DesiredGates_,DesiredStatesToSearch_,DesiredVertexNumber_]";
Begin["`Private`"];
GraphFinderExactVertex[InitialState_,DesiredGates_,DesiredStatesToSearch_,DesiredVertexNumber_]:=Block[{TotalGateSet = {},DesiredGateSet={},SingleQubits = Permutations[Range[Log[2,Length[InitialState]]],{1}],PairedQubits = Permutations[Range[Log[2,Length[InitialState]]],{2}],StateSet,SetOfStatesWithGates ,OldSet={},BooleanOff = True,StatesChecked ,CheckSet,TestState,RandomizedSequence,iii,Counter=0,jjj,lll,mmm},

(*Generate gate action desired.*)
For[jjj=1,jjj<Length[DesiredGates]+1,jjj++,
(*TotalGateSet = DeleteCases[TotalGateSet,DesiredGates[[kk]]];*)
If[DesiredGates[[jjj,1]]=="H",AppendTo[DesiredGateSet,HadamardAction[(Log[2,Length[InitialState]]),DesiredGates[[jjj,2,1]]]],If[DesiredGates[[jjj,1]]=="P",AppendTo[DesiredGateSet,PhaseAction[(Log[2,Length[InitialState]]),DesiredGates[[jjj,2,1]]]],AppendTo[DesiredGateSet,CNOTAction[(Log[2,Length[InitialState]]),DesiredGates[[jjj,2]]]],Print["Unable to interpret gate."]]]];

StateSet = NumberSelectedStateGenerator[InitialState,DesiredStatesToSearch];
RandomizedSequence = RandomSample[StateSet];

(*Iterate through all generated states.*)
For[iii=1,iii<Length[RandomizedSequence]+1,iii++,

(*Reset lists for upcoming While construct.*)
StatesChecked={};
OldSet={};
CheckSet = {RandomizedSequence[[iii]]};
SetOfStatesWithGates = {};

While[SubsetQ[StatesChecked,CheckSet]== False,

StatesChecked = DeleteDuplicates[Join[StatesChecked,CheckSet]];

For[lll=1,lll< Length[DesiredGateSet]+1,lll++,
For[mmm=1,mmm< Length[CheckSet]+1,mmm++,
TestState=DesiredGateSet[[lll]].CheckSet[[mmm]];
AppendTo[OldSet,((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState];AppendTo[SetOfStatesWithGates,{((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState,DesiredGates[[lll]],CheckSet[[mmm]],ReducedEntropyVectorBuilder[TestState]}]]];

CheckSet=DeleteDuplicates[OldSet];
OldSet={};
];
If[Length[StatesChecked]== DesiredVertexNumber,Print[RandomizedSequence[[iii]]];Print["There are: ",Length[DeleteDuplicates[ReducedEntropyVectors[StatesChecked]]]," different entropy vectors in this subgraph."];Print[SimpleReachabilityGraphBuilder[StatesChecked,DeleteDuplicates[SetOfStatesWithGates],DesiredGates]];Break[],If[Counter>10,Print["States Checked: ",iii,"/",Length[RandomizedSequence]];Counter=0,Counter++],Print["Confused."]];
]
];
End[];



GraphFinder::usage= "GraphFinder[InitialState_,DesiredGates_,DesiredStatesToSearch_,DesiredVertexNumber_]";
Begin["`Private`"];
GraphFinder[InitialState_,DesiredGates_,DesiredStatesToSearch_,DesiredVertexNumber_]:=Block[{TotalGateSet = {},DesiredGateSet={},SingleQubits = Permutations[Range[Log[2,Length[InitialState]]],{1}],PairedQubits = Permutations[Range[Log[2,Length[InitialState]]],{2}],StateSet,SetOfStatesWithGates ,OldSet={},BooleanOff = True,StatesChecked ,CheckSet,TestState,RandomizedSequence,iii,Counter=0,jjj,lll,mmm},

(*Generate gate action desired.*)
For[jjj=1,jjj<Length[DesiredGates]+1,jjj++,
(*TotalGateSet = DeleteCases[TotalGateSet,DesiredGates[[kk]]];*)
If[DesiredGates[[jjj,1]]=="H",AppendTo[DesiredGateSet,HadamardAction[(Log[2,Length[InitialState]]),DesiredGates[[jjj,2,1]]]],If[DesiredGates[[jjj,1]]=="P",AppendTo[DesiredGateSet,PhaseAction[(Log[2,Length[InitialState]]),DesiredGates[[jjj,2,1]]]],AppendTo[DesiredGateSet,CNOTAction[(Log[2,Length[InitialState]]),DesiredGates[[jjj,2]]]],Print["Unable to interpret gate."]]]];

StateSet = NumberSelectedStateGenerator[InitialState,DesiredStatesToSearch];
RandomizedSequence = RandomSample[StateSet];

(*Iterate through all generated states.*)
For[iii=1,iii<Length[RandomizedSequence]+1,iii++,

(*Reset lists for upcoming While construct.*)
StatesChecked={};
OldSet={};
CheckSet = {RandomizedSequence[[iii]]};
SetOfStatesWithGates = {};

While[SubsetQ[StatesChecked,CheckSet]== False,

StatesChecked = DeleteDuplicates[Join[StatesChecked,CheckSet]];

For[lll=1,lll< Length[DesiredGateSet]+1,lll++,
For[mmm=1,mmm< Length[CheckSet]+1,mmm++,
TestState=DesiredGateSet[[lll]].CheckSet[[mmm]];
AppendTo[OldSet,((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState];AppendTo[SetOfStatesWithGates,{((TestState[[FirstPosition[Flatten[Transpose[TestState]],_?(##!= 0&)][[1]]]][[1]])^(-1))*TestState,DesiredGates[[lll]],CheckSet[[mmm]],ReducedEntropyVectorBuilder[TestState]}]]];

CheckSet=DeleteDuplicates[OldSet];
OldSet={};
];
If[Length[StatesChecked]>DesiredVertexNumber,Print[RandomizedSequence[[iii]]];Print["There are: ",Length[DeleteDuplicates[EntropyVectors[StatesChecked]]]," in this subgraph."];Print[SimpleReachabilityGraphBuilder[StatesChecked,DeleteDuplicates[SetOfStatesWithGates],DesiredGates]];Break[],If[Counter>10,Print["States Checked: ",iii,"/",Length[RandomizedSequence]];Counter=0,Counter++],Print["Confused."]];
]
];
End[];



InteractiveReachabilityGraphBuilder::usage= "InteractiveReachabilityGraphBuilder[SetOfStabilizerStates_,SetOfStatesWithGates_,GatesToBeDisplayed_]";
Begin["`Private`"];
InteractiveReachabilityGraphBuilder[SetOfStabilizerStates_,SetOfStatesWithGates_,GatesToBeDisplayed_]:=Module[{NumberOfQubits = Log[2,Length[SetOfStabilizerStates[[1]]]],PossibleEntropyVectors = DeleteDuplicates[ReducedEntropyVectors[SetOfStabilizerStates]],GateBitPair={},HolographicStates = {},OperationIndices,EdgeSet={},VertexLabelSet={},VertexStyleSet = {},
VertexColorChoices = {Blue,Red,Orange,Purple,Yellow,Black,Cyan,Pink,Brown,LightGreen,White,LightRed,LightBlue,LightYellow,LightMagenta,LightOrange,LightGray,LightCyan,LightBrown},i,j,k,l,m,G},

If[Log[2,Length[SetOfStabilizerStates[[1]]]]>4,

For[m=1,m< Length[SetOfStabilizerStates]+1,m++,
AppendTo[VertexLabelSet,m-> Placed[Toggler[Invisible[m],{Style[ToString[Flatten[Transpose[SetOfStabilizerStates[[m]]]]],{Background-> LightYellow,CellFrameColor-> LightBlue,CellFrame-> 1}],Style[ToString[ReducedEntropyVectorBuilder[SetOfStabilizerStates[[m]]]],{Background-> LightBlue,CellFrameColor-> LightBlue,CellFrame-> 1}],Invisible[m]}],Center]]],

For[m=1,m< Length[SetOfStabilizerStates]+1,m++,
AppendTo[VertexLabelSet,m-> Placed[Toggler[Invisible[m],{Style[ToString[Flatten[Transpose[SetOfStabilizerStates[[m]]]]],{Background-> LightYellow,CellFrameColor-> LightBlue,CellFrame-> 1}],Style[ToString[ReducedEntropyVectorBuilder[SetOfStabilizerStates[[m]]]],{Background-> LightBlue,CellFrameColor-> LightBlue,CellFrame-> 1}],Invisible[m]}],Center]];AppendTo[VertexStyleSet ,m->VertexColorChoices[[Position[PossibleEntropyVectors,ReducedEntropyVectorBuilder[SetOfStabilizerStates[[m]]]][[1]]]][[1]]]]];

(*For the user's sanity*)
Print["Vertices colored and labeled."];

For[j=1,j<Length[SetOfStatesWithGates]+1,j++,AppendTo[GateBitPair,{ToString[SetOfStatesWithGates[[j,2]]],SetOfStatesWithGates[[j,3]]}]];
(*For the user's sanity*)
Print["Gates compiled."];

For[l=1,l< Length[GatesToBeDisplayed]+1,l++,
OperationIndices = Position[GateBitPair,GatesToBeDisplayed[[l]]];
For[k=1,k<Length[OperationIndices]+1,k++,
If[GatesToBeDisplayed[[l,1]]=="P" ,AppendTo[EdgeSet,DirectedEdge[Position[SetOfStabilizerStates,SetOfStatesWithGates[[OperationIndices[[k,1]],4]]][[1,1]],Position[SetOfStabilizerStates,SetOfStatesWithGates[[OperationIndices[[k,1]],1]]][[1,1]]]],If[GatesToBeDisplayed[[l,1]]=="H" ,AppendTo[EdgeSet,Style[Min[Position[SetOfStabilizerStates,SetOfStatesWithGates[[OperationIndices[[k,1]],4]]][[1,1]],Position[SetOfStabilizerStates,SetOfStatesWithGates[[OperationIndices[[k,1]],1]]][[1,1]]]<-> Max[Position[SetOfStabilizerStates,SetOfStatesWithGates[[OperationIndices[[k,1]],4]]][[1,1]],Position[SetOfStabilizerStates,SetOfStatesWithGates[[OperationIndices[[k,1]],1]]][[1,1]]],Blend[{Green,White},1-1/(SetOfStatesWithGates[[OperationIndices[[k,1]],3,1]])]]],AppendTo[EdgeSet,Style[Min[Position[SetOfStabilizerStates,SetOfStatesWithGates[[OperationIndices[[k,1]],4]]][[1,1]],Position[SetOfStabilizerStates,SetOfStatesWithGates[[OperationIndices[[k,1]],1]]][[1,1]]]<-> Max[Position[SetOfStabilizerStates,SetOfStatesWithGates[[OperationIndices[[k,1]],4]]][[1,1]],Position[SetOfStabilizerStates,SetOfStatesWithGates[[OperationIndices[[k,1]],1]]][[1,1]]],Blend[{Magenta,White},1-1/(SetOfStatesWithGates[[OperationIndices[[k,1]],3,1]])]]],Print["Unknown Operation."]]]];

(*For the user's sanity*)
Print["Edge set ",ToString[l]," of ",ToString[Length[GatesToBeDisplayed]], " added."]];

(*Uses the edge set compiled above to render the graph.*)
If[Log[2,Length[SetOfStabilizerStates[[1]]]]>4,G = Graph[DeleteDuplicates[EdgeSet],VertexLabels->VertexLabelSet],
G = Graph[DeleteDuplicates[EdgeSet],VertexLabels->VertexLabelSet,VertexStyle-> VertexStyleSet]];
G];
End[];



CycleFinder::usage= "CycleFinder[InitialGraph_,StartingBitAddress_,CycleLength_,ChosenCycleNumber_]";
Begin["`Private`"];
CycleFinder[InitialGraph_,StartingBitAddress_,CycleLength_,ChosenCycleNumber_]:=Module[{StateSet,ChosenCycle,CycleList,SortedCycleList,NumCycles,i,j},
If[Log[2,Length[StartingBitAddress]] == 1,StateSet = OneQubitStabilizerStates,If[Log[2,Length[StartingBitAddress]] == 2,StateSet = TwoQubitStabilizerStates,If[Log[2,Length[StartingBitAddress]] == 3,StateSet = ThreeQubitStabilizerStates,
If[Log[2,Length[StartingBitAddress]] == 4,StateSet = FourQubitStabilizerStates,
If[Log[2,Length[StartingBitAddress]] == 5,StateSet = FiveQubitStabilizerStates,Print["Invalid bit-address."],Print["Unable to interpret bit-address"]]]]]];
CycleList = FindCycle[{InitialGraph,Position[StateSet,Transpose[{StartingBitAddress}]][[1,1]]},{CycleLength},All];
NumCycles = Length[CycleList];
Print["There are ",NumCycles," cycles of length ", CycleLength, " involving ", StartingBitAddress, "." ];
For[i=1,i< Length[CycleList]+1,i++,
For[j=1,j<CycleLength+1,j++,CycleList[[i,j]]= Sort[CycleList[[i,j]]]]];
Print["Displaying cycle ",ChosenCycleNumber,":"];
HighlightGraph[InitialGraph,Style[CycleList[[ChosenCycleNumber]],Thick,Black]]];
End[];





(* ::Subsection:: *)
(*Subgraph Generators*)


ConnectedSubgraphGenerator::usage= "ConnectedSubgraphGenerator[InitialGraph_,DisplaySubgraphBoolean_]";
Begin["`Private`"];
ConnectedSubgraphGenerator[InitialGraph_,DisplaySubgraphBoolean_]:=Module[{UpdatedVertexList,SubgraphSet={},DropList={},VertexCountList={},VertexNumList={},SubgraphVertexList={},i,j,k,g,h},

UpdatedVertexList = VertexList[InitialGraph];
While[Length[UpdatedVertexList]>0,AppendTo[SubgraphSet,Flatten[ConnectedComponents[InitialGraph,UpdatedVertexList[[1]]]]];
UpdatedVertexList = DeleteCases[UpdatedVertexList,Alternatives@@Intersection[UpdatedVertexList,Flatten[ConnectedComponents[InitialGraph,UpdatedVertexList[[1]]]]]]];
For[i=1,i<Length[SubgraphSet]+1,i++,AppendTo[VertexCountList,Length[SubgraphSet[[i]]]]];
VertexNumList=DeleteDuplicates[VertexCountList];
Print["There are ",Length[SubgraphSet], " separate, connected subgraphs of vertex numbers ",VertexNumList ,"."];
For[j=1,j<Length[VertexNumList]+1,j++,Print["There are ",Count[VertexCountList,VertexNumList[[j]]]," connected subgraphs of ",VertexNumList[[j]]," vertices."]];

If[DisplaySubgraphBoolean== True,
For[g=1,g<Length[SubgraphSet]+1,g++,
For[h=1,h<Length[SubgraphSet]+1,h++,If[And[g!= h,IsomorphicGraphQ[Subgraph[InitialGraph,ConnectedComponents[InitialGraph,SubgraphSet[[g]]]],Subgraph[InitialGraph,ConnectedComponents[InitialGraph,SubgraphSet[[h]]]]]],AppendTo[DropList,{h}],Null,Print["Cannot determine if isomorphic."]
]];SubgraphSet = Delete[SubgraphSet,DropList];DropList={}];For[k=1,k<Length[SubgraphSet]+1,k++,Print[Subgraph[InitialGraph,ConnectedComponents[InitialGraph,SubgraphSet[[k]]]]]],Null]
];
End[];


SequenceGenerator::usage= "SequenceGenerator[InitialGraph_,GraphPath_,SetOfStates_,SetOfStatesWithEntropy_]";
Begin["`Private`"];
SequenceGenerator[InitialGraph_,GraphPath_,SetOfStates_,SetOfStatesWithEntropy_]:=Module[{FlattenedPath = Flatten[GraphPath],SequenceOfGates = {},SequenceOfQubits ={},SequenceOfOperations,PairsList={},i,j},
For[i=1,i<Length[FlattenedPath],i++,AppendTo[PairsList,First[Position[SetOfStatesWithEntropy,{SetOfStates[[FlattenedPath[[i+1]]]],_,_, SetOfStates[[FlattenedPath[[i]]]],_}]]]];
For[j=1,j<Length[Flatten[PairsList]]+1,j++,AppendTo[SequenceOfGates,SetOfStatesWithEntropy[[Flatten[PairsList][[j]]]][[2]]];AppendTo[SequenceOfQubits,SetOfStatesWithEntropy[[Flatten[PairsList][[j]]]][[3]]]];
SequenceOfOperations = {ReplaceAll[SequenceOfGates,{"H"-> HadamardGate,"CNOT"-> CNOTGate,"P"-> PhaseGate}],SequenceOfQubits};
SequenceOfOperations];
End[];



(* ::Subsection:: *)
(*Entropy Vectors*)


SubAdditivityCheck::usage= "SubAdditivityCheck[StateKet_]";
Begin["`Private`"];
SubAdditivityCheck[StateKet_]:=Module[{EV = EntropyVectorBuilder[StateKet],i},
For[i=1,i<Log[2,Length[StateKet]]+1,i++,If[,True,Print["Fails Subadditivity Check."]]]];
End[];


MutualInformationMeasure::usage= "MutualInformationMeasure[StateKet_,RegionA_,RegionB_]";
Begin["`Private`"];
MutualInformationMeasure[StateKet_,RegionA_,RegionB_]:= Module[{NormKet =Transpose[{Normalize[Flatten[Transpose[StateKet]]]}], MutualInformationMeasure},
MutualInformationMeasure = SvonNeumannBinary[NormKet,RegionA] + SvonNeumannBinary[NormKet,RegionB] - SvonNeumannBinary[NormKet,DeleteDuplicates[Join[RegionA,RegionB]]];
MutualInformationMeasure];
End[];


MMICheckerStateInput::usage= "MMICheckerStateInput[StateKet_]";
Begin["`Private`"];
MMICheckerStateInput[StateKet_]:=Module[{DisjointRegionTriplesList = {},RegionsList = Join[DeleteCases[Subsets[Range[Log[2,Length[StateKet]]],Log[2,Length[StateKet]]-1],{}],{Range[Log[2,Length[StateKet]]]}],RegionTriplesList,NumberOfFails = 0,i,j},

RegionTriplesList = Subsets[RegionsList,{3}];

For[i=1,i<Length[RegionTriplesList]+1,i++,If[And[DisjointQ[RegionTriplesList[[i,1]],RegionTriplesList[[i,2]]],DisjointQ[RegionTriplesList[[i,1]],RegionTriplesList[[i,3]]],DisjointQ[RegionTriplesList[[i,2]],RegionTriplesList[[i,3]]]],AppendTo[DisjointRegionTriplesList,RegionTriplesList[[i]]],Null]];

For[j=1,j<Length[DisjointRegionTriplesList]+1,j++,
If[MutualInformationMeasure[StateKet,DisjointRegionTriplesList[[j,1]],Union[RegionTriplesList[[j,2]],DisjointRegionTriplesList[[j,3]]]]>= MutualInformationMeasure[StateKet,DisjointRegionTriplesList[[j,1]],DisjointRegionTriplesList[[j,2]]]+MutualInformationMeasure[StateKet,DisjointRegionTriplesList[[j,1]],DisjointRegionTriplesList[[j,3]]],Null,NumberOfFails++;Print["Failure detected for: ",StateKet, " with entropy vector ", EntropyVectorBuilder[StateKet], " and regions ", DisjointRegionTriplesList[[j]]]]];

Print["Check Complete with ", NumberOfFails, " failures."]];
End[];



MMICheckerEntropyInput::usage= "MMICheckerEntropyInput[EntropyVector_]";
Begin["`Private`"];
MMICheckerEntropyInput[EntropyVector_]:=Module[{RegionSubsets = Join[DeleteCases[Subsets[Range[Log[2,Length[EntropyVector]+1]],Log[2,Length[EntropyVector]+1]-1],{}],{Range[Log[2,Length[EntropyVector]+1]]}],DisjointRegionTriplesList = {},RegionTriplesList = Subsets[DeleteCases[Subsets[Range[Log[2,Length[EntropyVector]+1]],3],{}],{3}],NumberOfFails = 0,i,j},

For[i=1,i<Length[RegionTriplesList]+1,i++,If[And[DisjointQ[RegionTriplesList[[i,1]],RegionTriplesList[[i,2]]],DisjointQ[RegionTriplesList[[i,1]],RegionTriplesList[[i,3]]],DisjointQ[RegionTriplesList[[i,2]],RegionTriplesList[[i,3]]]],AppendTo[DisjointRegionTriplesList,RegionTriplesList[[i]]],Null]];

For[j=1,j<Length[DisjointRegionTriplesList]+1,j++,
If[EntropyVector[[Position[RegionSubsets,Sort[Union[DisjointRegionTriplesList[[j,1]],DisjointRegionTriplesList[[j,2]]]]][[1,1]]]][[1]]+EntropyVector[[Position[RegionSubsets,Sort[Union[DisjointRegionTriplesList[[j,1]],DisjointRegionTriplesList[[j,3]]]]][[1,1]]]][[1]]+EntropyVector[[Position[RegionSubsets,Sort[Union[DisjointRegionTriplesList[[j,2]],DisjointRegionTriplesList[[j,3]]]]][[1,1]]]][[1]]>= EntropyVector[[Position[RegionSubsets,DisjointRegionTriplesList[[j,1]]][[1,1]]]][[1]] + EntropyVector[[Position[RegionSubsets,DisjointRegionTriplesList[[j,2]]][[1,1]]]][[1]] + EntropyVector[[Position[RegionSubsets,DisjointRegionTriplesList[[j,3]]][[1,1]]]] [[1]]+ EntropyVector[[Position[RegionSubsets,Sort[Union[DisjointRegionTriplesList[[j,1]],DisjointRegionTriplesList[[j,2]],DisjointRegionTriplesList[[j,3]]]]][[1,1]]]][[1]],Null,NumberOfFails++;Print["Failure detected for ",EntropyVector, " for regions ", DisjointRegionTriplesList[[j]]]]];

Print["Check Complete with ", NumberOfFails, " failures."]];
End[];


HolographicChecker::usage= "HolographicChecker[StateVector_]";
Begin["`Private`"];
HolographicChecker[StateVector_]:=Module[{SubsystemCount = Log[2,Length[Flatten[StateVector]]],EntropyVec = Flatten[EntropyVectorBuilder[StateVector]],Result=False},

(*There are no entangled single-party states, and all systems with sub-party number below 3 are holographic.*)
If[SubsystemCount<3,Result=True,

(*Verify the 3-party entropy inequalities (Strong Subadditivity and MMI) are satisfied.*)
If[SubsystemCount == 3,Result = And[EntropyVec[[4]]+EntropyVec[[6]]>= EntropyVec[[2]]+EntropyVec[[7]],EntropyVec[[4]]+EntropyVec[[5]]>= EntropyVec[[3]]+EntropyVec[[7]],EntropyVec[[5]]+EntropyVec[[6]]>= EntropyVec[[1]]+EntropyVec[[7]],EntropyVec[[4]]+EntropyVec[[5]]+EntropyVec[[6]]>= EntropyVec[[1]]+EntropyVec[[2]]+EntropyVec[[3]]+EntropyVec[[7]]],

(*Verify the 4-party entropy inequalities are satisfied.*)
If[SubsystemCount == 4,Result =False,Print["Subsystem count too high."]]]];
Result];
End[];


HolographicStateFinder::usage= "HolographicStateFinder[InitialGraph_,StabilizerStateSet_]";
Begin["`Private`"];
HolographicStateFinder[InitialGraph_,StabilizerStateSet_]:=Module[{HolographicVertices={},i},For[i=1,i<Length[VertexList[InitialGraph]]+1,i++,If[HolographicChecker[StabilizerStateSet[[i]]]== True,AppendTo[HolographicVertices,i],Null]];
HighlightGraph[InitialGraph,Style[HolographicVertices,Large,Green]]];
End[];


(* ::Section:: *)
(**)
