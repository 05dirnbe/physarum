/*
* Author: Michael Dirnberger <mtd@mpi-inf.mpg.de>
*
* Copyright (c)  Max-Planck-Institute Saarbruecken (Germany)
*
* This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
* WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*
*/
/**
* @file
* @author  Michael Dirnberger <mtd@mpi-inf.mpg.de>
* @version 1.0
* @ingroup
*
* Some description should go here, cmon, dont be lazy
*
**/

#include <iostream>

template < typename BigFloat = double >
class LinearSolver_Default {

 public:

    LinearSolver_Default() {  std::cout << "Constructing LinearSolver_Default" << std::endl; }
    ~LinearSolver_Default() {  std::cout<<"Destructing LinearSolver_Default"<<std::endl; }
//  ~LinearSolver_Default() { }

  protected:
     void solve() const { std::cout<<"\tDefault linear solver here"<<std::endl;}

};

template < typename BigFloat = double >
class LinearSolver_A {

 public:

  LinearSolver_A() {  std::cout<<"Constructing LinearSolver_A"<<std::endl; }
  ~LinearSolver_A() {  std::cout<<"Destructing LinearSolver_A"<<std::endl; }
//  ~LinearSolver_A() { }

  protected:
     void solve() const { std::cout<<"\tLinear solver A here"<<std::endl;}

};





template < typename BigFloat = long double >
class MidPointRule {

 public:

  MidPointRule() {  std::cout<<"Constructing MidPointRule"<<std::endl; }
  ~MidPointRule() {  std::cout<<"Destructing MidPointRule"<<std::endl; }
//  ~MidPointRule() { }

  protected:
     void differentiate() const { std::cout<<"\tUsing Midpoint rule"<<std::endl;}

};

template < typename BigFloat = long double >
class ThreePointRule {

 public:

  ThreePointRule() {  std::cout<<"Constructing ThreePointRule"<<std::endl; }
  ~ThreePointRule() {  std::cout<<"Destructing ThreePointRule"<<std::endl; }
//  ~ThreePointRule() { }

  protected:
     void differentiate() const { std::cout<<"\tUsing ThreePointRule rule"<<std::endl;}

};



template
<
        typename BigFloat = long double,
        class NumericDifferentiationPolicy = MidPointRule< BigFloat >
>
class PhysarumUpdate_Default : public NumericDifferentiationPolicy {

    typedef NumericDifferentiationPolicy Differentiation;

    public:

        PhysarumUpdate_Default() {  std::cout<<"Constructing PhysarumUpdate_Default"<<std::endl; }
        ~PhysarumUpdate_Default() {  std::cout<<"Destructing PhysarumUpdate_Default"<<std::endl; }
//      ~PhysarumUpdate_Default() { }


    private:
        //here goes stuff only this class itself needs to know about

    protected:

        void update() const {

            std::cout<<"\tDefault Physarum Update here"<<std::endl;
            std::cout<<"\t";
            Differentiation::differentiate();
        }

};

template
<
        class BigFloat = long double,
        class NumericDifferentiationPolicy = MidPointRule< BigFloat >
>
class PhysarumUpdate_Special : public NumericDifferentiationPolicy {

    typedef NumericDifferentiationPolicy Differentiation;

    public:

        PhysarumUpdate_Special() {  std::cout<<"Constructing PhysarumUpdate_Special"<<std::endl; }
        ~PhysarumUpdate_Special() {  std::cout<<"Destructing PhysarumUpdate_Special"<<std::endl; }
//      ~PhysarumUpdate_Special() { }


    private:
        //here goes stuff only this class itself needs to know about

    protected:

        void update() const {

            std::cout<<"\tSpecial Physarum Update here"<<std::endl;
            std::cout<<"\t";
            Differentiation::differentiate();
        }

};







//here goes the solver
template
<
            class BigFloat = long double,
            class PhysarumUpdatePolicy = PhysarumUpdate_Default< BigFloat >,
            class LinearSolverPolicy = LinearSolver_Default< BigFloat >
>
class PhysarumSolver : public PhysarumUpdatePolicy, LinearSolverPolicy {

    typedef PhysarumUpdatePolicy UpdateRule;
    typedef LinearSolverPolicy LinearSolver;

    public:


        PhysarumSolver() {  std::cout<<"Constructing Physarum solver"<<std::endl; }
        ~PhysarumSolver() {  std::cout<<"Destructing Physarum solver"<<std::endl; }
//        ~PhysarumSolver() {  }

        void run() const { compute_potentials(); compute_flux(); }

    private:

        void compute_potentials() const { LinearSolver::solve(); }        //solve a linear system
        void compute_flux() const { UpdateRule::update(); }              //apply the physarum specific stuff

};


using namespace std;

int main() {


    //note: non-empty base class destructor will cause failure of inlining as it stands
    //note: not sure if base class requires virtual destructors

    PhysarumSolver< long double, PhysarumUpdate_Special< ThreePointRule<double> > , LinearSolver_A<float> > Solver;
//    PhysarumSolver<> Solver;

    Solver.run();
}
