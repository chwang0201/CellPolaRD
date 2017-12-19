#include "Grid/grid_dist_id.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "Vector/vector_dist.hpp"
#include "data_type/aggregate.hpp"
#include "Grid/Iterators/grid_dist_id_iterator_dec_skin.hpp"
#include "timer.hpp"
#include <math.h>

//! \cond [constants] \endcond

constexpr int U = 0;
constexpr int V = 1;
constexpr int W = 2;
constexpr int Z = 3;
constexpr int N = 4;
constexpr int O = 5;
constexpr int R = 6;
constexpr int F = 7;

constexpr int x = 0;
constexpr int y = 1;

//! \cond [constants] \endcond

//! \cond [init fun] \endcond

void init(grid_dist_id<2,double,aggregate<double,double,double,double,double,double,double,double> > & Old, grid_dist_id<2,double,aggregate<double,double,double,double,double,double,double,double> > & New, Box<2,double> & domain)
{

//! \cond [init fun] \endcond

	//! \cond [init uv] \endcond

	double UM = 0.25; //Concentration of Free E-Cad #/grid
    double VM = 0.25; //Concentration of Free target #/grid
    double WM = 0.01;
    double ZM = 0.025;
    double NM = 0.25;
    double OM = 0.25;
    double RM = 0.01;

	auto it = Old.getDomainIterator();

    /*double noiseUM [100*100] = {};
    double noiseVM [100*100] = {};
    double noiseZM [100*100] = {};
    double noiseNM [100*100] = {};
    double noiseOM [100*100] = {};
    double noiseRM [100*100] = {};

    for (int j = 0; j < 100*100; j++)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> d(UM,0.01);
            double zeta = d(gen);
            noiseUM [j] = zeta; 
        }

    for (int j = 0; j < 100*100; j++)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> d(NM,0.01);
            double zeta = d(gen);
            noiseNM [j] = zeta; 
        }

    for (int j = 0; j < 100*100; j++)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> d(VM,0.01);
            double zeta = d(gen);
            noiseVM [j] = zeta; 
        }

    for (int j = 0; j < 100*100; j++)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> d(ZM,0.001);
            double zeta = d(gen);
            noiseZM [j] = zeta; 
        }

    for (int j = 0; j < 100*100; j++)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> d(OM,0.0001);
            double zeta = d(gen);
            noiseOM [j] = zeta; 
        }

    for (int j = 0; j < 100*100; j++)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> d(RM,0.0001);
            double zeta = d(gen);
            noiseRM [j] = zeta; 
        }*/

	while (it.isNext())
	{
		// Get the local grid key
		auto key = it.get();
        auto key_g = Old.getGKey(key);

        /*int axis1 = floor ((((double)std::rand())/RAND_MAX)*100*100);
        int axis2 = floor ((((double)std::rand())/RAND_MAX)*100*100);
        int axis3 = floor ((((double)std::rand())/RAND_MAX)*100*100);
        int axis4 = floor ((((double)std::rand())/RAND_MAX)*100*100);
        int axis5 = floor ((((double)std::rand())/RAND_MAX)*100*100);
        int axis6 = floor ((((double)std::rand())/RAND_MAX)*100*100);*/

        Old.template get<U>(key) = UM + UM*((((double)std::rand())/RAND_MAX - 0.5)/20.0); //#/um^2 
        Old.template get<V>(key) = VM + VM*((((double)std::rand())/RAND_MAX - 0.5)/20.0); //#/um^2 
        Old.template get<W>(key) = WM + WM*((((double)std::rand())/RAND_MAX - 0.5)/20.0); //#/um^2 
        Old.template get<Z>(key) = ZM; //+ ZM*((((double)std::rand())/RAND_MAX - 0.5));
        //ZM - (0.0075 - (0.0075/2500)*(key_g.get(0)-0.5*100.0)*(key_g.get(0)-0.5*100.0) + 0.0075 - (0.0075/2500)*(key_g.get(1)-0.5*100.0)*(key_g.get(1)-0.5*100.0));
        
        Old.template get<N>(key) = NM + NM*((((double)std::rand())/RAND_MAX - 0.5)/20.0); //#/um^2 
        Old.template get<O>(key) = OM + OM*((((double)std::rand())/RAND_MAX - 0.5)/20.0); //#/um^2 
        Old.template get<R>(key) = RM + RM*((((double)std::rand())/RAND_MAX - 0.5)/20.0); //#/um^2 
        Old.template get<F>(key) = Old.template get<Z>(key);

		New.template get<U>(key) = 0.0;
		New.template get<V>(key) = 0.0;
		New.template get<W>(key) = 0.0;
		New.template get<Z>(key) = 0.0;
        New.template get<N>(key) = 0.0;
        New.template get<O>(key) = 0.0;
        New.template get<R>(key) = 0.0;
		New.template get<F>(key) = New.template get<F>(key);

		++it;
	}

	grid_key_dx<2> start({(long int)std::floor(Old.size(0)*2.00f/domain.getHigh(0)),(long int)std::floor(Old.size(1)*2.00f/domain.getHigh(1))});
    grid_key_dx<2> stop ({(long int)std::ceil (Old.size(0)*8.00f/domain.getHigh(0)),(long int)std::ceil (Old.size(1)*8.00f/domain.getHigh(1))});
    auto it_init = Old.getSubDomainIterator(start,stop);
    
    while (it_init.isNext())
    {
        auto key = it_init.get();
        auto key_g = Old.getGKey(key);
        
        Old.template get<Z>(key) = ZM - (0.0075 - (0.0075/900)*(key_g.get(0)-0.5*100.0)*(key_g.get(0)-0.5*100.0) + 0.0075 - (0.0075/900)*(key_g.get(1)-0.5*100.0)*(key_g.get(1)-0.5*100.0));

        //Old.template get<Z>(key) = ZM + 7.5 - sqrt (7.5*7.5 - (key_g.get(0)-0.5*100.0)*(key_g.get(0)-0.5*100.0) + (key_g.get(1)-0.5*100.0)*(key_g.get(1)-0.5*100.0));
        //Old.template get<Z>(key) = ZM + ZM*((((double)std::rand())/RAND_MAX - 0.5)/50.0);

        //std::cout << key_g.get(0) << std::endl;
        
        ++it_init;
    }

}

/*Define boundary condition. Create a imaginary boundary at the actual boundary. And mirror the value at the actual boundary to have no flux flowing out*/
void mirror2(grid_dist_id<2, double, aggregate<double,double,double,double,double,double,double,double,double>> & g_dist)
{
    Box<2,size_t> A({0,0},{g_dist.size(0)-1,g_dist.size(1)-1});
    Box<2,size_t> B = A;
    
    size_t bc[2] = {NON_PERIODIC,NON_PERIODIC};
    
    grid_dist_id_iterator_dec_skin<CartDecomposition<2,double>> it_dec(g_dist.getDecomposition(),g_dist.getGridInfo(),A,B,bc);
    
    
    while (it_dec.isNext() == true)
    {
        auto key = it_dec.get_int();
        auto p = g_dist.getGKey(key);
        
        if (p.get(x) == 0)
        {
            if (g_dist.get<Z>(key) - 0.027 > 0)
            {
                g_dist.get<Z>(key.move(x,-1)) = g_dist.get<Z>(key) - 0.0015*0.02;
                g_dist.get<Z>(key.move(x,-2)) = g_dist.get<Z>(key.move(x,-1));
                g_dist.get<U>(key.move(x,-1)) = g_dist.get<U>(key);
                g_dist.get<U>(key.move(x,-2)) = g_dist.get<U>(key.move(x,-1));
                g_dist.get<V>(key.move(x,-1)) = g_dist.get<V>(key);
                g_dist.get<V>(key.move(x,-2)) = g_dist.get<V>(key.move(x,-1));
                g_dist.get<W>(key.move(x,-1)) = g_dist.get<W>(key);
                g_dist.get<W>(key.move(x,-2)) = g_dist.get<W>(key.move(x,-1));
                g_dist.get<N>(key.move(x,-1)) = g_dist.get<N>(key);
                g_dist.get<N>(key.move(x,-2)) = g_dist.get<N>(key.move(x,-1));
                g_dist.get<O>(key.move(x,-1)) = g_dist.get<O>(key);
                g_dist.get<O>(key.move(x,-2)) = g_dist.get<O>(key.move(x,-1));
                g_dist.get<R>(key.move(x,-1)) = g_dist.get<R>(key);
                g_dist.get<R>(key.move(x,-2)) = g_dist.get<R>(key.move(x,-1));
            }
            else if (g_dist.get<Z>(key) - 0.027 <= 0)
            {
                g_dist.get<Z>(key.move(x,-1)) = g_dist.get<Z>(key);
                g_dist.get<Z>(key.move(x,-2)) = g_dist.get<Z>(key.move(x,-1));
                g_dist.get<U>(key.move(x,-1)) = g_dist.get<U>(key);
                g_dist.get<U>(key.move(x,-2)) = g_dist.get<U>(key.move(x,-1));
                g_dist.get<V>(key.move(x,-1)) = g_dist.get<V>(key);
                g_dist.get<V>(key.move(x,-2)) = g_dist.get<V>(key.move(x,-1));
                g_dist.get<W>(key.move(x,-1)) = g_dist.get<W>(key);
                g_dist.get<W>(key.move(x,-2)) = g_dist.get<W>(key.move(x,-1));
                g_dist.get<N>(key.move(x,-1)) = g_dist.get<N>(key);
                g_dist.get<N>(key.move(x,-2)) = g_dist.get<N>(key.move(x,-1));
                g_dist.get<O>(key.move(x,-1)) = g_dist.get<O>(key);
                g_dist.get<O>(key.move(x,-2)) = g_dist.get<O>(key.move(x,-1));
                g_dist.get<R>(key.move(x,-1)) = g_dist.get<R>(key);
                g_dist.get<R>(key.move(x,-2)) = g_dist.get<R>(key.move(x,-1));
            }
            
        }
        else if (p.get(x) == g_dist.size(x) - 1)
        {
            if (g_dist.get<Z>(key) - 0.027 > 0)
            {
                g_dist.get<Z>(key.move(x,1)) = g_dist.get<Z>(key) - 0.0015*0.02;
                g_dist.get<Z>(key.move(x,2)) = g_dist.get<Z>(key.move(x,1));
                g_dist.get<U>(key.move(x,1)) = g_dist.get<U>(key);
                g_dist.get<U>(key.move(x,2)) = g_dist.get<U>(key.move(x,1));
                g_dist.get<V>(key.move(x,1)) = g_dist.get<V>(key);
                g_dist.get<V>(key.move(x,2)) = g_dist.get<V>(key.move(x,1));
                g_dist.get<W>(key.move(x,1)) = g_dist.get<W>(key);
                g_dist.get<W>(key.move(x,2)) = g_dist.get<W>(key.move(x,1));
                g_dist.get<N>(key.move(x,1)) = g_dist.get<N>(key);
                g_dist.get<N>(key.move(x,2)) = g_dist.get<N>(key.move(x,1));
                g_dist.get<O>(key.move(x,1)) = g_dist.get<O>(key);
                g_dist.get<O>(key.move(x,2)) = g_dist.get<O>(key.move(x,1));
                g_dist.get<R>(key.move(x,1)) = g_dist.get<R>(key);
                g_dist.get<R>(key.move(x,2)) = g_dist.get<R>(key.move(x,1));
            }
            else if (g_dist.get<Z>(key) - 0.027 <= 0)
            {
                g_dist.get<Z>(key.move(x,1)) = g_dist.get<Z>(key);
                g_dist.get<Z>(key.move(x,2)) = g_dist.get<Z>(key.move(x,1));
                g_dist.get<U>(key.move(x,1)) = g_dist.get<U>(key);
                g_dist.get<U>(key.move(x,2)) = g_dist.get<U>(key.move(x,1));
                g_dist.get<V>(key.move(x,1)) = g_dist.get<V>(key);
                g_dist.get<V>(key.move(x,2)) = g_dist.get<V>(key.move(x,1));
                g_dist.get<W>(key.move(x,1)) = g_dist.get<W>(key);
                g_dist.get<W>(key.move(x,2)) = g_dist.get<W>(key.move(x,1));
                g_dist.get<N>(key.move(x,1)) = g_dist.get<N>(key);
                g_dist.get<N>(key.move(x,2)) = g_dist.get<N>(key.move(x,1));
                g_dist.get<O>(key.move(x,1)) = g_dist.get<O>(key);
                g_dist.get<O>(key.move(x,2)) = g_dist.get<O>(key.move(x,1));
                g_dist.get<R>(key.move(x,1)) = g_dist.get<R>(key);
                g_dist.get<R>(key.move(x,2)) = g_dist.get<R>(key.move(x,1));
            }
            
        }
        if (p.get(y) == 0)
        {
            if (g_dist.get<Z>(key) - 0.027 > 0)
            {
                g_dist.get<Z>(key.move(y,-1)) = g_dist.get<Z>(key) - 0.0015*0.02;
                g_dist.get<Z>(key.move(y,-2)) = g_dist.get<Z>(key.move(y,-1));  
                g_dist.get<U>(key.move(y,-1)) = g_dist.get<U>(key);
                g_dist.get<U>(key.move(y,-2)) = g_dist.get<U>(key.move(y,-1)); 
                g_dist.get<V>(key.move(y,-1)) = g_dist.get<V>(key);
                g_dist.get<V>(key.move(y,-2)) = g_dist.get<V>(key.move(y,-1)); 
                g_dist.get<W>(key.move(y,-1)) = g_dist.get<W>(key);
                g_dist.get<W>(key.move(y,-2)) = g_dist.get<W>(key.move(y,-1)); 
                g_dist.get<N>(key.move(y,-1)) = g_dist.get<N>(key);
                g_dist.get<N>(key.move(y,-2)) = g_dist.get<N>(key.move(y,-1)); 
                g_dist.get<O>(key.move(y,-1)) = g_dist.get<O>(key);
                g_dist.get<O>(key.move(y,-2)) = g_dist.get<O>(key.move(y,-1)); 
                g_dist.get<R>(key.move(y,-1)) = g_dist.get<R>(key);  
                g_dist.get<R>(key.move(y,-2)) = g_dist.get<R>(key.move(y,-1)); 
            }
            else if (g_dist.get<Z>(key) - 0.027 <= 0)
            {
                g_dist.get<Z>(key.move(y,-1)) = g_dist.get<Z>(key);
                g_dist.get<Z>(key.move(y,-2)) = g_dist.get<Z>(key.move(y,-1));  
                g_dist.get<U>(key.move(y,-1)) = g_dist.get<U>(key);
                g_dist.get<U>(key.move(y,-2)) = g_dist.get<U>(key.move(y,-1)); 
                g_dist.get<V>(key.move(y,-1)) = g_dist.get<V>(key);
                g_dist.get<V>(key.move(y,-2)) = g_dist.get<V>(key.move(y,-1)); 
                g_dist.get<W>(key.move(y,-1)) = g_dist.get<W>(key);
                g_dist.get<W>(key.move(y,-2)) = g_dist.get<W>(key.move(y,-1)); 
                g_dist.get<N>(key.move(y,-1)) = g_dist.get<N>(key);
                g_dist.get<N>(key.move(y,-2)) = g_dist.get<N>(key.move(y,-1)); 
                g_dist.get<O>(key.move(y,-1)) = g_dist.get<O>(key);
                g_dist.get<O>(key.move(y,-2)) = g_dist.get<O>(key.move(y,-1)); 
                g_dist.get<R>(key.move(y,-1)) = g_dist.get<R>(key);  
                g_dist.get<R>(key.move(y,-2)) = g_dist.get<R>(key.move(y,-1)); 
            }
            
        }
        else if (p.get(y) == g_dist.size(y) - 1)
        {
            if (g_dist.get<Z>(key) - 0.027 > 0)
            {
                g_dist.get<Z>(key.move(y,1)) = g_dist.get<Z>(key) - 0.0015*0.02;
                g_dist.get<Z>(key.move(y,2)) = g_dist.get<Z>(key.move(y,1));  
                g_dist.get<U>(key.move(y,1)) = g_dist.get<U>(key);
                g_dist.get<U>(key.move(y,2)) = g_dist.get<U>(key.move(y,1));
                g_dist.get<V>(key.move(y,1)) = g_dist.get<V>(key);
                g_dist.get<V>(key.move(y,2)) = g_dist.get<V>(key.move(y,1));
                g_dist.get<W>(key.move(y,1)) = g_dist.get<W>(key);
                g_dist.get<W>(key.move(y,2)) = g_dist.get<W>(key.move(y,1));
                g_dist.get<N>(key.move(y,1)) = g_dist.get<N>(key);
                g_dist.get<N>(key.move(y,2)) = g_dist.get<N>(key.move(y,1));
                g_dist.get<O>(key.move(y,1)) = g_dist.get<O>(key);
                g_dist.get<O>(key.move(y,2)) = g_dist.get<O>(key.move(y,1));
                g_dist.get<R>(key.move(y,1)) = g_dist.get<R>(key); 
                g_dist.get<R>(key.move(y,2)) = g_dist.get<R>(key.move(y,1));
            }
            else if (g_dist.get<Z>(key) - 0.027 <= 0)
            {
                g_dist.get<Z>(key.move(y,1)) = g_dist.get<Z>(key);
                g_dist.get<Z>(key.move(y,2)) = g_dist.get<Z>(key.move(y,1));  
                g_dist.get<U>(key.move(y,1)) = g_dist.get<U>(key);
                g_dist.get<U>(key.move(y,2)) = g_dist.get<U>(key.move(y,1));
                g_dist.get<V>(key.move(y,1)) = g_dist.get<V>(key);
                g_dist.get<V>(key.move(y,2)) = g_dist.get<V>(key.move(y,1));
                g_dist.get<W>(key.move(y,1)) = g_dist.get<W>(key);
                g_dist.get<W>(key.move(y,2)) = g_dist.get<W>(key.move(y,1));
                g_dist.get<N>(key.move(y,1)) = g_dist.get<N>(key);
                g_dist.get<N>(key.move(y,2)) = g_dist.get<N>(key.move(y,1));
                g_dist.get<O>(key.move(y,1)) = g_dist.get<O>(key);
                g_dist.get<O>(key.move(y,2)) = g_dist.get<O>(key.move(y,1));
                g_dist.get<R>(key.move(y,1)) = g_dist.get<R>(key); 
                g_dist.get<R>(key.move(y,2)) = g_dist.get<R>(key.move(y,1));
            }
           
        } 
        if (p.get(x) == 0 && p.get(y) == 0)
        {
            g_dist.get<Z>(key.move(x,-1).move(y,-1)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-2).move(y,-2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-2).move(y,-1)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-1).move(y,-2)) = g_dist.get<Z>(key);
            g_dist.get<U>(key.move(x,-1).move(y,-1)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-2).move(y,-2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-2).move(y,-1)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-1).move(y,-2)) = g_dist.get<U>(key);
            g_dist.get<V>(key.move(x,-1).move(y,-1)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-2).move(y,-2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-2).move(y,-1)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-1).move(y,-2)) = g_dist.get<V>(key);
            g_dist.get<W>(key.move(x,-1).move(y,-1)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-2).move(y,-2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-2).move(y,-1)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-1).move(y,-2)) = g_dist.get<W>(key);
            g_dist.get<N>(key.move(x,-1).move(y,-1)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-2).move(y,-2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-2).move(y,-1)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-1).move(y,-2)) = g_dist.get<N>(key);
            g_dist.get<O>(key.move(x,-1).move(y,-1)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-2).move(y,-2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-2).move(y,-1)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-1).move(y,-2)) = g_dist.get<O>(key);
            g_dist.get<R>(key.move(x,-1).move(y,-1)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-2).move(y,-2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-2).move(y,-1)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-1).move(y,-2)) = g_dist.get<R>(key);
        }
        else if (p.get(x) == 0 && p.get(y) == g_dist.size(y) - 1)
        {
            g_dist.get<Z>(key.move(x,-1).move(y,1)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-2).move(y,2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-1).move(y,2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-2).move(y,1)) = g_dist.get<Z>(key);
            g_dist.get<U>(key.move(x,-1).move(y,1)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-2).move(y,2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-1).move(y,2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-2).move(y,1)) = g_dist.get<U>(key);
            g_dist.get<V>(key.move(x,-1).move(y,1)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-2).move(y,2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-1).move(y,2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-2).move(y,1)) = g_dist.get<V>(key);
            g_dist.get<W>(key.move(x,-1).move(y,1)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-2).move(y,2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-1).move(y,2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-2).move(y,1)) = g_dist.get<W>(key);
            g_dist.get<N>(key.move(x,-1).move(y,1)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-2).move(y,2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-1).move(y,2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-2).move(y,1)) = g_dist.get<N>(key);
            g_dist.get<O>(key.move(x,-1).move(y,1)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-2).move(y,2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-1).move(y,2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-2).move(y,1)) = g_dist.get<O>(key);
            g_dist.get<R>(key.move(x,-1).move(y,1)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-2).move(y,2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-1).move(y,2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-2).move(y,1)) = g_dist.get<R>(key);
        }
        else if (p.get(x) == g_dist.size(x) - 1 && p.get(y) == 0)
        {
            g_dist.get<Z>(key.move(x,1).move(y,-1)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,2).move(y,-2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,1).move(y,-2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,2).move(y,-1)) = g_dist.get<Z>(key);
            g_dist.get<U>(key.move(x,1).move(y,-1)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,2).move(y,-2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,1).move(y,-2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,2).move(y,-1)) = g_dist.get<U>(key);
            g_dist.get<V>(key.move(x,1).move(y,-1)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,2).move(y,-2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,1).move(y,-2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,2).move(y,-1)) = g_dist.get<V>(key);
            g_dist.get<W>(key.move(x,1).move(y,-1)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,2).move(y,-2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,1).move(y,-2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,2).move(y,-1)) = g_dist.get<W>(key);
            g_dist.get<N>(key.move(x,1).move(y,-1)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,2).move(y,-2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,1).move(y,-2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,2).move(y,-1)) = g_dist.get<N>(key);
            g_dist.get<O>(key.move(x,1).move(y,-1)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,2).move(y,-2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,1).move(y,-2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,2).move(y,-1)) = g_dist.get<O>(key);
            g_dist.get<R>(key.move(x,1).move(y,-1)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,2).move(y,-2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,1).move(y,-2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,2).move(y,-1)) = g_dist.get<R>(key);
        }
        else if (p.get(x) == g_dist.size(x) - 1 && p.get(y) == g_dist.size(y) - 1)
        {
            g_dist.get<Z>(key.move(x,1).move(y,1)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,2).move(y,2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,1).move(y,2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,2).move(y,1)) = g_dist.get<Z>(key);
            g_dist.get<U>(key.move(x,1).move(y,1)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,2).move(y,2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,1).move(y,2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,2).move(y,1)) = g_dist.get<U>(key);
            g_dist.get<V>(key.move(x,1).move(y,1)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,2).move(y,2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,1).move(y,2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,2).move(y,1)) = g_dist.get<V>(key);
            g_dist.get<W>(key.move(x,1).move(y,1)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,2).move(y,2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,1).move(y,2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,2).move(y,1)) = g_dist.get<W>(key);
            g_dist.get<N>(key.move(x,1).move(y,1)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,2).move(y,2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,1).move(y,2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,2).move(y,1)) = g_dist.get<N>(key);
            g_dist.get<O>(key.move(x,1).move(y,1)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,2).move(y,2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,1).move(y,2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,2).move(y,1)) = g_dist.get<O>(key);
            g_dist.get<R>(key.move(x,1).move(y,1)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,2).move(y,2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,1).move(y,2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,2).move(y,1)) = g_dist.get<R>(key);
        }
        ++it_dec;
    }
}

void mirror(grid_dist_id<2, double, aggregate<double,double,double,double,double,double,double,double>> & g_dist)
{
    Box<2,size_t> A({0,0},{g_dist.size(0)-1,g_dist.size(1)-1});
    Box<2,size_t> B = A;
    
    size_t bc[2] = {NON_PERIODIC,NON_PERIODIC};
    
    grid_dist_id_iterator_dec_skin<CartDecomposition<2,double>> it_dec(g_dist.getDecomposition(),g_dist.getGridInfo(),A,B,bc);
    
    
    while (it_dec.isNext() == true)
    {
        auto key = it_dec.get_int();
        auto p = g_dist.getGKey(key);

           if (p.get(x) == 0)
            {
            		//g_dist.get<Z>(key) = g_dist.get<F>(key);
                    g_dist.get<Z>(key.move(x,-1)) = g_dist.get<Z>(key);
                    g_dist.get<Z>(key.move(x,-2)) = g_dist.get<Z>(key.move(x,-1));
                    g_dist.get<U>(key.move(x,-1)) = g_dist.get<U>(key);
                    g_dist.get<U>(key.move(x,-2)) = g_dist.get<U>(key.move(x,-1));
                    g_dist.get<V>(key.move(x,-1)) = g_dist.get<V>(key);
                    g_dist.get<V>(key.move(x,-2)) = g_dist.get<V>(key.move(x,-1));
                    g_dist.get<W>(key.move(x,-1)) = g_dist.get<W>(key);
                    g_dist.get<W>(key.move(x,-2)) = g_dist.get<W>(key.move(x,-1));
                    g_dist.get<N>(key.move(x,-1)) = g_dist.get<N>(key);
                    g_dist.get<N>(key.move(x,-2)) = g_dist.get<N>(key.move(x,-1));
                    g_dist.get<O>(key.move(x,-1)) = g_dist.get<O>(key);
                    g_dist.get<O>(key.move(x,-2)) = g_dist.get<O>(key.move(x,-1));
                    g_dist.get<R>(key.move(x,-1)) = g_dist.get<R>(key);
                    g_dist.get<R>(key.move(x,-2)) = g_dist.get<R>(key.move(x,-1));
            
            }
            else if (p.get(x) == g_dist.size(x) - 1)
            {

            		//g_dist.get<Z>(key) = g_dist.get<F>(key);
                    g_dist.get<Z>(key.move(x,1)) = g_dist.get<Z>(key);
                    g_dist.get<Z>(key.move(x,2)) = g_dist.get<Z>(key.move(x,1));
                    g_dist.get<U>(key.move(x,1)) = g_dist.get<U>(key);
                    g_dist.get<U>(key.move(x,2)) = g_dist.get<U>(key.move(x,1));
                    g_dist.get<V>(key.move(x,1)) = g_dist.get<V>(key);
                    g_dist.get<V>(key.move(x,2)) = g_dist.get<V>(key.move(x,1));
                    g_dist.get<W>(key.move(x,1)) = g_dist.get<W>(key);
                    g_dist.get<W>(key.move(x,2)) = g_dist.get<W>(key.move(x,1));
                    g_dist.get<N>(key.move(x,1)) = g_dist.get<N>(key);
                    g_dist.get<N>(key.move(x,2)) = g_dist.get<N>(key.move(x,1));
                    g_dist.get<O>(key.move(x,1)) = g_dist.get<O>(key);
                    g_dist.get<O>(key.move(x,2)) = g_dist.get<O>(key.move(x,1));
                    g_dist.get<R>(key.move(x,1)) = g_dist.get<R>(key);
                    g_dist.get<R>(key.move(x,2)) = g_dist.get<R>(key.move(x,1));
           
            }
            if (p.get(y) == 0)
            {

            		//g_dist.get<Z>(key) = g_dist.get<F>(key);
                    g_dist.get<Z>(key.move(y,-1)) = g_dist.get<Z>(key);
                    g_dist.get<Z>(key.move(y,-2)) = g_dist.get<Z>(key.move(y,-1));  
                    g_dist.get<U>(key.move(y,-1)) = g_dist.get<U>(key);
                    g_dist.get<U>(key.move(y,-2)) = g_dist.get<U>(key.move(y,-1)); 
                    g_dist.get<V>(key.move(y,-1)) = g_dist.get<V>(key);
                    g_dist.get<V>(key.move(y,-2)) = g_dist.get<V>(key.move(y,-1)); 
                    g_dist.get<W>(key.move(y,-1)) = g_dist.get<W>(key);
                    g_dist.get<W>(key.move(y,-2)) = g_dist.get<W>(key.move(y,-1)); 
                    g_dist.get<N>(key.move(y,-1)) = g_dist.get<N>(key);
                    g_dist.get<N>(key.move(y,-2)) = g_dist.get<N>(key.move(y,-1)); 
                    g_dist.get<O>(key.move(y,-1)) = g_dist.get<O>(key);
                    g_dist.get<O>(key.move(y,-2)) = g_dist.get<O>(key.move(y,-1)); 
                    g_dist.get<R>(key.move(y,-1)) = g_dist.get<R>(key);  
                    g_dist.get<R>(key.move(y,-2)) = g_dist.get<R>(key.move(y,-1)); 

            }
            else if (p.get(y) == g_dist.size(y) - 1)
            {
            		//g_dist.get<Z>(key) = g_dist.get<F>(key);
                    g_dist.get<Z>(key.move(y,1)) = g_dist.get<Z>(key);
                    g_dist.get<Z>(key.move(y,2)) = g_dist.get<Z>(key.move(y,1));  
                    g_dist.get<U>(key.move(y,1)) = g_dist.get<U>(key);
                    g_dist.get<U>(key.move(y,2)) = g_dist.get<U>(key.move(y,1));
                    g_dist.get<V>(key.move(y,1)) = g_dist.get<V>(key);
                    g_dist.get<V>(key.move(y,2)) = g_dist.get<V>(key.move(y,1));
                    g_dist.get<W>(key.move(y,1)) = g_dist.get<W>(key);
                    g_dist.get<W>(key.move(y,2)) = g_dist.get<W>(key.move(y,1));
                    g_dist.get<N>(key.move(y,1)) = g_dist.get<N>(key);
                    g_dist.get<N>(key.move(y,2)) = g_dist.get<N>(key.move(y,1));
                    g_dist.get<O>(key.move(y,1)) = g_dist.get<O>(key);
                    g_dist.get<O>(key.move(y,2)) = g_dist.get<O>(key.move(y,1));
                    g_dist.get<R>(key.move(y,1)) = g_dist.get<R>(key); 
                    g_dist.get<R>(key.move(y,2)) = g_dist.get<R>(key.move(y,1));
            } 
            if (p.get(x) == 0 && p.get(y) == 0)
            {
            g_dist.get<Z>(key.move(x,-1).move(y,-1)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-2).move(y,-2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-2).move(y,-1)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-1).move(y,-2)) = g_dist.get<Z>(key);
            g_dist.get<U>(key.move(x,-1).move(y,-1)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-2).move(y,-2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-2).move(y,-1)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-1).move(y,-2)) = g_dist.get<U>(key);
            g_dist.get<V>(key.move(x,-1).move(y,-1)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-2).move(y,-2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-2).move(y,-1)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-1).move(y,-2)) = g_dist.get<V>(key);
            g_dist.get<W>(key.move(x,-1).move(y,-1)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-2).move(y,-2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-2).move(y,-1)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-1).move(y,-2)) = g_dist.get<W>(key);
            g_dist.get<N>(key.move(x,-1).move(y,-1)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-2).move(y,-2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-2).move(y,-1)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-1).move(y,-2)) = g_dist.get<N>(key);
            g_dist.get<O>(key.move(x,-1).move(y,-1)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-2).move(y,-2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-2).move(y,-1)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-1).move(y,-2)) = g_dist.get<O>(key);
            g_dist.get<R>(key.move(x,-1).move(y,-1)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-2).move(y,-2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-2).move(y,-1)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-1).move(y,-2)) = g_dist.get<R>(key);
            }
            else if (p.get(x) == 0 && p.get(y) == g_dist.size(y) - 1)
            {
            g_dist.get<Z>(key.move(x,-1).move(y,1)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-2).move(y,2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-1).move(y,2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,-2).move(y,1)) = g_dist.get<Z>(key);
            g_dist.get<U>(key.move(x,-1).move(y,1)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-2).move(y,2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-1).move(y,2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,-2).move(y,1)) = g_dist.get<U>(key);
            g_dist.get<V>(key.move(x,-1).move(y,1)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-2).move(y,2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-1).move(y,2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,-2).move(y,1)) = g_dist.get<V>(key);
            g_dist.get<W>(key.move(x,-1).move(y,1)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-2).move(y,2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-1).move(y,2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,-2).move(y,1)) = g_dist.get<W>(key);
            g_dist.get<N>(key.move(x,-1).move(y,1)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-2).move(y,2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-1).move(y,2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,-2).move(y,1)) = g_dist.get<N>(key);
            g_dist.get<O>(key.move(x,-1).move(y,1)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-2).move(y,2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-1).move(y,2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,-2).move(y,1)) = g_dist.get<O>(key);
            g_dist.get<R>(key.move(x,-1).move(y,1)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-2).move(y,2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-1).move(y,2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,-2).move(y,1)) = g_dist.get<R>(key);
            }
            else if (p.get(x) == g_dist.size(x) - 1 && p.get(y) == 0)
            {
            g_dist.get<Z>(key.move(x,1).move(y,-1)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,2).move(y,-2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,1).move(y,-2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,2).move(y,-1)) = g_dist.get<Z>(key);
            g_dist.get<U>(key.move(x,1).move(y,-1)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,2).move(y,-2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,1).move(y,-2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,2).move(y,-1)) = g_dist.get<U>(key);
            g_dist.get<V>(key.move(x,1).move(y,-1)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,2).move(y,-2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,1).move(y,-2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,2).move(y,-1)) = g_dist.get<V>(key);
            g_dist.get<W>(key.move(x,1).move(y,-1)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,2).move(y,-2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,1).move(y,-2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,2).move(y,-1)) = g_dist.get<W>(key);
            g_dist.get<N>(key.move(x,1).move(y,-1)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,2).move(y,-2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,1).move(y,-2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,2).move(y,-1)) = g_dist.get<N>(key);
            g_dist.get<O>(key.move(x,1).move(y,-1)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,2).move(y,-2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,1).move(y,-2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,2).move(y,-1)) = g_dist.get<O>(key);
            g_dist.get<R>(key.move(x,1).move(y,-1)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,2).move(y,-2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,1).move(y,-2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,2).move(y,-1)) = g_dist.get<R>(key);
            }
            else if (p.get(x) == g_dist.size(x) - 1 && p.get(y) == g_dist.size(y) - 1)
            {
            g_dist.get<Z>(key.move(x,1).move(y,1)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,2).move(y,2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,1).move(y,2)) = g_dist.get<Z>(key);
            g_dist.get<Z>(key.move(x,2).move(y,1)) = g_dist.get<Z>(key);
            g_dist.get<U>(key.move(x,1).move(y,1)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,2).move(y,2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,1).move(y,2)) = g_dist.get<U>(key);
            g_dist.get<U>(key.move(x,2).move(y,1)) = g_dist.get<U>(key);
            g_dist.get<V>(key.move(x,1).move(y,1)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,2).move(y,2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,1).move(y,2)) = g_dist.get<V>(key);
            g_dist.get<V>(key.move(x,2).move(y,1)) = g_dist.get<V>(key);
            g_dist.get<W>(key.move(x,1).move(y,1)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,2).move(y,2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,1).move(y,2)) = g_dist.get<W>(key);
            g_dist.get<W>(key.move(x,2).move(y,1)) = g_dist.get<W>(key);
            g_dist.get<N>(key.move(x,1).move(y,1)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,2).move(y,2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,1).move(y,2)) = g_dist.get<N>(key);
            g_dist.get<N>(key.move(x,2).move(y,1)) = g_dist.get<N>(key);
            g_dist.get<O>(key.move(x,1).move(y,1)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,2).move(y,2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,1).move(y,2)) = g_dist.get<O>(key);
            g_dist.get<O>(key.move(x,2).move(y,1)) = g_dist.get<O>(key);
            g_dist.get<R>(key.move(x,1).move(y,1)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,2).move(y,2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,1).move(y,2)) = g_dist.get<R>(key);
            g_dist.get<R>(key.move(x,2).move(y,1)) = g_dist.get<R>(key);
            }         
        
        ++it_dec;
    }
}
//! \cond [end fun] \endcond


int main(int argc, char* argv[])
{

	//! \cond [init lib] \endcond

	openfpm_init(&argc,&argv);

	// domain
	Box<2,double> domain({0.0,0.0},{10.0,10.0});
	
	// grid size
	size_t sz[2] = {100,100};
    int x_axis = sz[0];
    int y_axis = sz[1];

    double spacing_ = domain.getLow(0) / (sz[0] - 1);
    int GridNum = x_axis*y_axis;

	// Define periodicity of the grid
	periodicity<2> bc = {NON_PERIODIC,NON_PERIODIC};
	
    if (bc.bc[0] == NON_PERIODIC)
    {
        spacing_ = domain.getHigh(0) / (sz[0] - 1);
    }
    else
    {
        spacing_ = domain.getHigh(0) / (sz[0]);
    }

	// Ghost in grid unit
	//Ghost<2,double> g(10/(sz[0]-10));
	Ghost<2,double> g(spacing_ * 2.000001);

	// Diffusion constant for specie 
	double du = 0.5;
	double dv = 0.5;
	double dw = du*1e-5;
    double dn = 0.5;
    double dO = dn;
    double dr = dn*1e-5;

	//double Kd  = 2.125*1e-7;  // 1#/um2
	double Kon1 = 300;
    double Koff1 = 1;    //Expand or decrease effect of unbinding, 1/s
    double Kon4 = 300;
    double Koff4 = 1; 
    double lambda1 = 10;     //Bending energy of bounded junction
    double lambda2 = 100;
    double kappa = 200;      //bending regidity
    double gamma = 3.1;      //surface tension
    double M = 1e-5;
    double Zw1 = 0.02;
    double Zw2 = 0.01;          //Length of bounded pair um
    double delTa1 = 0.0025;
    double delTa2 = 0.005;
    double Kb = 1.38*1e-23;
    double T = 310;
    double gammaj = gamma/(Kb*T);
    double Kbond1 = 40*1e-6;
    double Kbond2 = 40*1e-6;
    double Kon2 = 0.03;
    double Kon3 = 0.006;
    double Koff3 = 0.0012;
    double tho2 = 0.1;
    double thooff2 = 0.01;

	grid_dist_id<2, double, aggregate<double,double,double,double,double,double,double,double>> Old(sz,domain,g,bc);

	// New grid with the decomposition of the old grid
	grid_dist_id<2, double, aggregate<double,double,double,double,double,double,double,double>> New(Old.getDecomposition(),sz,domain,g);

	// spacing of the grid on x and y
	double spacing[2] = {Old.spacing(0),Old.spacing(1)};

    // Number of timesteps
	double deltaT = 
	0.5*(spacing[x]*spacing[x]*spacing[x]*spacing[x])/(8*M*kappa);
	//0.5*(spacing[x]*spacing[x])/(4*du);

    double Time = 7200;

	size_t timeSteps = Time/deltaT;

    int OutPutTime = floor (1/deltaT);

	std::cout << timeSteps << std::endl;
    std::cout << deltaT << std::endl;

	//! \cond [init grid] \endcond

	//! \cond [init uvc] \endcond

	init(Old,New,domain);

    Old.ghost_get<U,V,W,Z,N,O,R,F>();
    mirror(Old);

	/*Old.write("Initial");
    
    mirror(Old);
    
    Old.write("After");
 
    openfpm_finalize();
    
    return 0;*/

	//Old.write("test");

	size_t count = 0;
	//size_t counter = 0;
	Old.template ghost_get<U,V,W,Z,N,O,R,F>();

	// because we assume that spacing[x] == spacing[y] we use formula 2
	// and we calculate the prefactor of Eq 2

    double noise [100*100] = {};

    for (int j = 0; j < GridNum; j++)
        {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> d(0,0.01);
            double zeta = d(gen);
            noise [j] = zeta; 
        }

	for (size_t i = 0; i < timeSteps; ++i)
	{
		auto it = Old.getDomainIterator();


        //std::clock_t start;
        //double duration;

        //start = std::clock();

		//counter = 0;

		while (it.isNext())
		{

            auto key = it.get();
            auto k = Old.getGKey(key);

            double dist = 0.01 * sqrt((k.get(x) - 50)*(k.get(x) - 50) + (k.get(y) - 50)*(k.get(y) - 50));

            double DiffuWeigh = (exp(-10*dist));
            double DiffuWeigh2 = tanh(3*abs(dist - 0.7071));
            double dud = du;//*DiffuWeigh;
            double dvd = dv;//*DiffuWeigh;
            double dwd = dw;//*DiffuWeigh;
            double dnd = dn;//*DiffuWeigh;
            double dOd = dO;//*DiffuWeigh;
            double drd = dr;

            double uFactor = deltaT * dud/(spacing[x]*spacing[x]);
            double vFactor = deltaT * dvd/(spacing[x]*spacing[x]);
            double wFactor = deltaT * dwd/(spacing[x]*spacing[x]);
            double nFactor = deltaT * dnd/(spacing[x]*spacing[x]);
            double oFactor = deltaT * dOd/(spacing[x]*spacing[x]);
            double rFactor = deltaT * drd/(spacing[x]*spacing[x]);

            double KonZ1 = Kon1*(exp (-(Old.get<Z>(key) - Zw1)*(Old.get<Z>(key) - Zw1)/(2*delTa1*delTa1)));
            double KoffZ1 = Koff1*(exp ((Old.get<Z>(key) - Zw1)*(Old.get<Z>(key) - Zw1)/(2*delTa2*delTa2)));
            double KonZ2 = Kon4*(exp (-(Old.get<Z>(key) - Zw2)*(Old.get<Z>(key) - Zw2)/(2*delTa1*delTa1))); 
            double KoffZ2 = Koff4*(exp ((Old.get<Z>(key) - Zw2)*(Old.get<Z>(key) - Zw2)/(2*delTa2*delTa2))); 
            double Kon2Z = Kon2*(exp (-(Old.get<O>(key))));
            double Koff3Z = Koff3*(1 - exp ((Old.get<O>(key))));
            double Kon3Z = Kon3*(exp (-(Old.get<O>(key)))); 
            double thooffZ2 = thooff2*(exp ((Old.get<Z>(key) - Zw2)*(Old.get<Z>(key) - Zw2)/(2*delTa2*delTa2)));

            double Fw = lambda1*(Old.get<Z>(key) - Zw1);
            double Gw = (Old.get<W>(key.move(x,1)) - Old.get<W>(key))*(Old.get<Z>(key.move(x,1)) - Old.get<Z>(key)) +
                        (Old.get<W>(key.move(y,1)) - Old.get<W>(key))*(Old.get<Z>(key.move(y,1)) - Old.get<Z>(key));
            double Mw = Old.get<W>(key)*(Old.get<Z>(key.move(x,1)) + Old.get<Z>(key.move(x,-1)) + Old.get<Z>(key.move(y,1)) + Old.get<Z>(key.move(y,-1)) - 4*Old.get<Z>(key));

            double Fr = lambda2*(Old.get<Z>(key) - Zw2);
            double Gr = (Old.get<R>(key.move(x,1)) - Old.get<R>(key))*(Old.get<Z>(key.move(x,1)) - Old.get<Z>(key)) +
                        (Old.get<R>(key.move(y,1)) - Old.get<R>(key))*(Old.get<Z>(key.move(y,1)) - Old.get<Z>(key));
            double Mr = Old.get<R>(key)*(Old.get<Z>(key.move(x,1)) + Old.get<Z>(key.move(x,-1)) + Old.get<Z>(key.move(y,1)) + Old.get<Z>(key.move(y,-1)) - 4*Old.get<Z>(key));

            double FourthOrder = (Old.get<Z>(key.move(x,2)) - 4.0*Old.get<Z>(key.move(x,1)) + 6.0*Old.get<Z>(key) - 4.0*Old.get<Z>(key.move(x,-1)) + Old.get<Z>(key.move(x,-2)) + 
            					  Old.get<Z>(key.move(y,2)) - 4.0*Old.get<Z>(key.move(y,1)) + 6.0*Old.get<Z>(key) - 4.0*Old.get<Z>(key.move(y,-1)) + Old.get<Z>(key.move(y,-2)))/(spacing[x]*spacing[x]*spacing[x]*spacing[x]);
            double SecondOrder = (Old.get<Z>(key.move(x,1)) + Old.get<Z>(key.move(x,-1)) + Old.get<Z>(key.move(y,1)) + Old.get<Z>(key.move(y,-1)) - 4.0*Old.get<Z>(key))/(spacing[x]*spacing[x]);
            int axis = floor ((((double)std::rand())/RAND_MAX)*GridNum);

            double RUVW = KonZ1 * Old.get<U>(key) * Old.get<V>(key) - KoffZ1 * Old.get<W>(key);
            double RWNO = Kon2 * Old.get<W>(key) * Old.get<N>(key) * Old.get<N>(key) + Kon3 * Old.get<W>(key) * Old.get<N>(key) - Koff3 * Old.get<O>(key);
            double RNOR = KonZ2 * Old.get<N>(key) * Old.get<O>(key) - KoffZ2 * Old.get<R>(key);

            //double KonZ1 = Kon*(exp (-Kbond1*(Old.get<Z>(key) - Zw1)*(Old.get<Z>(key) - Zw1)/(2*Kb*T)));
            //double KoffZ1 = Koff*(exp (Kbond1*(Old.get<Z>(key) - Zw1)*(Old.get<Z>(key) - Zw1)/(2*Kb*T)));
            //double KonZ2 = KLon*(exp (-Kbond2*(Old.get<Z>(key) - Zw2)*(Old.get<Z>(key) - Zw2)/(2*Kb*T)));
            //double KoffZ2 = KLoff*(exp (Kbond2*(Old.get<Z>(key) - Zw2)*(Old.get<Z>(key) - Zw2)/(2*Kb*T)));

            //std::cout << counter << "		" << FourthOrder << std::endl;
            //std::cout << axis << std::endl;
            //std::cout << std::setprecision (15) << noise[axis] << std::endl;

            /*std::random_device rd;
            std::mt19937 gen(rd());
            std::normal_distribution<> d(0,0.1);
            double zeta = d(gen);*/

			New.get<U>(key) = Old.get<U>(key) + uFactor * (
										Old.get<U>(key.move(x,1)) +
										Old.get<U>(key.move(x,-1)) +
										Old.get<U>(key.move(y,1)) +
										Old.get<U>(key.move(y,-1)) -
										4.0*Old.get<U>(key)) +
                                        deltaT * (-RUVW);

			New.get<V>(key) = Old.get<V>(key) + vFactor * (
										Old.get<V>(key.move(x,1)) +
										Old.get<V>(key.move(x,-1)) +
										Old.get<V>(key.move(y,1)) +
										Old.get<V>(key.move(y,-1)) -
										4*Old.get<V>(key)) +
                                        deltaT * (-RUVW);
  			//New.get<test>(key) = 6;
  			New.get<W>(key) = Old.get<W>(key) + wFactor * (
                                        Old.get<W>(key.move(x,1)) +
                                        Old.get<W>(key.move(x,-1)) +
                                        Old.get<W>(key.move(y,1)) +
                                        Old.get<W>(key.move(y,-1)) -
                                        4.0*Old.get<W>(key) + Fw * (Gw + Mw)) +
                                        deltaT * (RUVW);

            New.get<N>(key) = Old.get<N>(key) + nFactor * (
                                        Old.get<N>(key.move(x,1)) +
                                        Old.get<N>(key.move(x,-1)) +
                                        Old.get<N>(key.move(y,1)) +
                                        Old.get<N>(key.move(y,-1)) -
                                        4.0*Old.get<N>(key)) +
                                        deltaT * (-RNOR);

            // update based on Eq 2
            New.get<O>(key) = Old.get<O>(key) + oFactor * (
                                        Old.get<O>(key.move(x,1)) +
                                        Old.get<O>(key.move(x,-1)) +
                                        Old.get<O>(key.move(y,1)) +
                                        Old.get<O>(key.move(y,-1)) -
                                        4*Old.get<O>(key)) + 
                                        deltaT * (-RNOR);

            New.get<R>(key) = Old.get<R>(key) + rFactor * (
                                        Old.get<R>(key.move(x,1)) +
                                        Old.get<R>(key.move(x,-1)) +
                                        Old.get<R>(key.move(y,1)) +
                                        Old.get<R>(key.move(y,-1)) -
                                        4.0*Old.get<R>(key) + Fr * (Gr + Mr)) +
                                        deltaT * (RNOR);
 
            
            New.get<Z>(key) = Old.get<Z>(key) - deltaT * M *(lambda1 * Old.get<W>(key) * (Old.get<Z>(key) - Zw1) +
                                                             lambda2 * Old.get<R>(key) * (Old.get<Z>(key) - Zw2) -
                                                            gamma*SecondOrder + kappa*FourthOrder) /*+ deltaT*noise[axis]*/;

            New.get<F>(key) = Old.get<F>(key);
            //New.get<F>(key) = noise[axis];
            //std::cout << New.get<Z>(key) << std::endl;
			++it;
			//counter = counter + 1;
		}

        if (i % (OutPutTime) == 0)
        {
            Old.write("Output_15_ParaMem_0GradBC_Slow",count);
            std::cout << i << std::endl;
            count++;
        }

		Old.copy(New);

		Old.ghost_get<U,V,W,Z,N,O,R,F>();
        mirror(Old);


        //duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

        //std::cout<<"printf: "<< duration <<'\n';

	}
	openfpm_finalize();
}
