import java.lang.Math;
import java.util.*;

class GeneticAlgo {
    
    public static double[][] newPop = new double[200][20];
    public static double[] fittest = new double[20];
    public static double[] secFittest = new double[20];
    public static double[] newSolu = new double[20];
    public static double fit;
    public static double prevFit;
    public static double mutationRate = 0.0001;
    public static Random rn = new Random();
    public static int count = 0;
    
    public static void main(String[] args) {
        // Start time of the algorithm when ran.
        long startT = System.currentTimeMillis();

        // initialize a new population.
        double[][] population = createPop();
        // get the fittest solution from population.
        fit = Assess.getTest1(getFittestSol(population));

        // while the fittest solution is above the value 0 then loop.
        while (fit > 0) {
            // Checks to see if the prevfit and current fittest values are fairly similar.
            if ((prevFit - fit) < 0.00000000001 && fit > 1) {
                population = createPop();
            }
            // Checks to see if the prevfit and fit are the same value.
            if (prevFit == fit && fit < 1) {
                population = mutateAll(population);
            }

            count += 1;

            // do a selection of the two best results.
            selection(population);

            prevFit = fit;

            // calculate new best fitness
            fit = Assess.getTest1(getFittestSol(population));

            System.out.println("Generation " + count + " Fittest : " + fit);
        }
        
        // End time of algorithm when completed.
        long endT = System.currentTimeMillis();
        // Prints the time taken for the genetic algorithm to execute.
        System.out.println("Total execution time was: " + ((endT - startT) / 1000.0) + " seconds");

    }

    //////////////////// Genetic Algorithm ///////////////////

    // Creates a new population.
    public static double[][] createPop() {
        double[][] pop = new double[200][20];
        for (int i = 0; i < 200; i++) {
            pop[i] = createSol();
        }
        return pop;
    }

    // Creates a solution.
    public static double[] createSol() {
        double[] sol = new double[20];
        for (int i = 0; i < sol.length; i++) {
            sol[i] = Math.random() * Math.round(5.12 * (Math.random() - Math.random()));
        }
        return sol;
    }

    // Provides a selection process for the Genetic algorithm.
    public static double[][] selection(double[][] pop) {
        // Selects the fittest and second fittest solutions from the population
        fittest = getFittestSol(pop);
        secFittest = getSecondFittestSol(pop);
        // Iterates through the current population and performs mutations and
        // crossovers.
        for (int i = 0; i < pop.length / 2; i++) {
            // on the random change rand is true then the fittest solution is added to the
            // new population, else second fittest is added.
            boolean rand = (Math.random() > 0.5);
            if (rand) {
                pop[i] = mutate(fittest);
            } else {
                pop[i] = mutate(secFittest);
            }

            // perform a crossover on the two fittest solutions in the current population.
            double[] newSol = crossover(fittest, secFittest);
            // mutate the cross over.
            mutate(newSol);
            // add mutated crossover to the new population.
            pop[199 - i] = newSol;
        }
        return pop;
    }

    // finds the best fitness.
    public static double[] getFittestSol(double[][] pop) {
        double[] fittest = new double[20];
        double fit = Assess.getTest1(pop[0]);
        for (int i = 0; i < pop.length; i++) {
            if (Assess.getTest1(pop[i]) <= fit) {
                fit = Assess.getTest1(pop[i]);
                fittest = pop[i];
            }
        }
        return fittest;
    }

    // finds the second best solution.
    public static double[] getSecondFittestSol(double[][] pop) {
        double[] secFittest = new double[20];
        double fit1 = Assess.getTest1(pop[0]);
        double fit2 = Assess.getTest1(pop[0]);
        for (int i = 0; i < pop.length; i++) {
            if (Assess.getTest1(pop[i]) <= fit1) {
                fit1 = Assess.getTest1(pop[i]);
            }
            if (Assess.getTest1(pop[i]) <= fit2 && Assess.getTest1(pop[i]) > fit1) {
                fit2 = Assess.getTest1(pop[i]);
                secFittest = pop[i];
            }
        }
        return secFittest;
    }

    // Mutates all values in the population.
    public static double[][] mutateAll(double[][] pop) {
        for (int i = 0; i < pop.length - 1; i++) {
            pop[i] = mutate(pop[i]);
        }
        return pop;
    }

    // Mutates a selected candidate solution.
    public static double[] mutate(double[] sol) {
        // checks to see if the fittest value of the new population has not changed
        // significantly.
        if ((prevFit - fit) < 0.0000001 && fit < 1) {
            // change the mutation rate to a more precise number
            mutationRate = 0.0000001;
            // checks the fittest solution needs a more precise mutation rate.
            if (prevFit - fit < 0.000000001) {
                mutationRate = 0.0000000001;
            }
        }

        // loops through the solution and mutates the selected data.
        for (int i = 0; i < sol.length; i++) {
            boolean sel = new Random().nextBoolean();
            // A random integer selected between the values -5 and 5
            int num = new Random().nextInt(5 - (-5)) + (-5);
            // The mutation that will be added or subtracted from the solution data
            double mutation = num * mutationRate;
            // on the random chance sel is true then add mutation, if not then subtract.
            if (sel == true) {
                // Checks to see if the new mutated data is within the range of -5 and 5.
                if ((sol[i] + mutation) < 5 && (sol[i] + mutation) > -5) {
                    sol[i] = sol[i] + mutation;
                } else {
                    sol[i] = sol[i];
                }
            } else {
                if ((sol[i] - mutation) < 5 && (sol[i] - mutation) > -5) {
                    sol[i] = sol[i] - mutation;
                } else {
                    sol[i] = sol[i];
                }
            }
        }
        return sol;
    }

    // Performs a cross over on two arrays.
    public static double[] crossover(double[] sol1, double[] sol2) {
        double[] cross = new double[20];
        double[] tail = new double[10];
        double[] head = new double[10];
        Random rnd = new Random();
        // Loops through the two solutions and adds data at random indexes to tail and
        // head arrays.
        for (int i = 0; i < 10; i++) {
            tail[i] = sol1[rnd.nextInt(20)];
            head[i] = sol2[rnd.nextInt(20)];
        }

        // Combines tail and head arrays together into the final cross array.
        System.arraycopy(tail, 0, cross, 0, 10);
        System.arraycopy(head, 0, cross, 10, 10);

        return cross;
    }
}
