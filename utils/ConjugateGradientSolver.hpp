/*
 * This source file is part of FDTD UNIUD.
 *
 * Copyright (C) 2017, Bernard Kapidani - kapidani.bernard@spes.uniud.it
 * Department of Electrical Engineering, University of Udine
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Udine nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR(s) ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(s) BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
 
 //file Simulation.hpp
 #include "Utilities.hpp"

class ConjugateGradientSolver
{
    Eigen::SparseMatrix<double>  A_;
   Eigen::VectorXd pre_;
   double tol, final_error;
   uint32_t N_, iter, maxiter;

   public:
   ConjugateGradientSolver()
   {
      tol = 1e-8;
      maxiter = 5000;
   }

    /* Conjugate gradient, as in
     * "Templates for the Solution of Linear Systems:
     *  Building Blocks for Iterative Methods"
     */
    Eigen::VectorXd solveWithGuess(const Eigen::VectorXd& b, const Eigen::VectorXd& resultx)
    {
        assert(b.size() == N_);
      assert(resultx.size() == N_);
      
        double beta;
        double rhoA_ = 0, rho_p;
        
        this->iter = 0;
        double residual0, residual;
        double a;
        Eigen::VectorXd r(N_), z(N_), p(N_), q(N_);
      auto x = resultx;
      
      timecounter tc;
      tc.tic();
        r = b - A_*x;
        tc.toc();
      
      // std::cout << "Multiplying by mass matrix takes " << tc << " seconds" << std::endl;
      
      residual0 = residual = r.norm();
        while ( residual/residual0 > tol && iter < maxiter )
        {
            z = pre_.cwiseProduct(r);
            //z = r;
            rho_p = rhoA_;
            rhoA_ = r.dot(z);
            
            if ( iter == 0 )
                p = z;
            else
            {
                beta = rhoA_/rho_p;
                p = z + beta * p;
            }
            
            q = A_*p;
            
            a = rhoA_/(p.dot(q));
            
            x = x + a*p;
            r = r - a*q;
            
            residual = r.norm();
            
            iter++;
        }
      
        this->final_error = residual/residual0;
        return x;
    }
   
   uint32_t iterations(void) const { return iter; }
   double   error(void) const { return final_error; }
   void setMaxIterations(const uint32_t& maxit) { this->maxiter= maxit; }
   void setTolerance(const double& t) { this->tol = t; }

    void compute(Eigen::SparseMatrix<double>& mb)
    {
      assert(mb.rows() == mb.cols());
      
      mb.makeCompressed();
      this->A_ = std::move(mb);
        this->N_ = A_.rows();
        this->pre_ = Eigen::VectorXd::Zero(N_);
      
      for(uint32_t j=0; j<A_.outerSize(); ++j)
      {
         Eigen::SparseMatrix<double>::InnerIterator it(A_,j);
         while(it && it.index() !=j)
            ++it;
         if (it && it.index()==j && it.value() != double(0))
            pre_(j) = double(1)/it.value();
         else
            pre_(j) = double(1);
      }
    }

};
