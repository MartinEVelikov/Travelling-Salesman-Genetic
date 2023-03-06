#include <iostream>
#include <cmath>
#include <limits>
#include <chrono>
#include <vector>
#include <string>
#include <random>
#include <algorithm>

using namespace std;

const int FLOAT_MIN = 0;
const int FLOAT_MAX = 1;
int generateRandomInteger(int size) {
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(1, (size - 1));
    return distrib(gen);
}

double generateRandomDouble() {
    std::random_device rd;
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<float> distr(FLOAT_MIN, FLOAT_MAX);
    return distr(eng);
}
bool areSame(double a, double b) {
    return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}

class City
{
private:
    double x;
    double y;
public:

    City();
    City(double x, double y);
    double getX() const;
    double getY() const;
    double distanceTo(const City& city) const;
    bool operator==(const City& rhs);
    friend ostream& operator<<(ostream& os, City& city)
    {
        os << "( " << city.getX() << ", " << city.getY() << ")";
        return os;
    }
};
City::City()
{
    this->x = generateRandomDouble() * 200;
    this->y = generateRandomDouble() * 200;
}

City::City(double x, double y)
{
    this->x = x;
    this->y = y;
}

double City::getX() const
{
    return this->x;
}

double City::getY() const
{
    return this->y;
}

double City::distanceTo(const City& city) const
{
    double xDistance = fabs(getX() - city.getX());
    double yDistance = fabs(getY() - city.getY());
    double distance = sqrt((xDistance * xDistance) + (yDistance * yDistance));

    return distance;
}

bool City::operator==(const City& rhs) {
    if (areSame(this->x, rhs.x) && areSame(this->y, rhs.y)) return true;
    else return false;
}

static vector<City>destinationCities;
const double mutationRate = 0.015;
const int tournamentSize = 5;

class Tour
{
private:
    vector<City> tour;
    double fitness = 0.0;
    double distance = 0.0;
public:
    Tour();
    Tour(vector<City> tour);
    void generateIndividual();
    City getCity(const int& tourPosition) const;
    void setCity(int tourPosition, City city);
    double getFitness();
    double getDistance();
    int tourSize() const;
    bool containsCity(const City& city);
    friend ostream& operator<<(ostream& os, const Tour& tour);
};
Tour::Tour()
{
    for (int i = 0; i < destinationCities.size(); i++)
    {
        tour.push_back(City(0.0, 0.0));
    }
}

Tour::Tour(vector<City> tour)
{
    this->tour = tour;
}


void Tour::generateIndividual()
{
    for (int cityIndex = 0; cityIndex < destinationCities.size(); cityIndex++)
    {
        setCity(cityIndex, destinationCities[cityIndex]);
    }
    random_shuffle(tour.begin(), tour.end());

}


City Tour::getCity(const int& tourPosition) const
{
    return tour[tourPosition];
}


void Tour::setCity(int tourPosition, City city)
{
    tour[tourPosition] = city;
    fitness = 0.0;
    distance = 0.0;
}


double Tour::getFitness() 
{
    if (fitness == 0.0)
    {
        fitness = 1.0 / getDistance();
    }
    return fitness;
}


double Tour::getDistance()
{
    double tourDistance = 0.0;
    for (int cityIndex = 0; cityIndex < tourSize() - 1; cityIndex++)
    {
        City fromCity = getCity(cityIndex);
        City destinationCity;
        destinationCity = getCity(cityIndex + 1);
        tourDistance += fromCity.distanceTo(destinationCity);
    }
    distance = tourDistance;
    return distance;
}


int Tour::tourSize() const
{
    return tour.size();
}


bool Tour::containsCity(const City& city)
{
    if (std::count(tour.begin(), tour.end(), city))
    {
        return true;
    }
    else
    {
        return false;
    }
}

ostream& operator<<(ostream& os, const Tour& tour)
{
    cout << '|';
    for (auto city : tour.tour)
    {
        cout << city << '|';
    }
    return os;
}

class Population
{
private:
    vector<Tour> tours;
public:
    Population(int populationSize, bool initialise);
    void saveTour(int index,const Tour& tour);
    Tour getTour(int index) const;
    Tour getFittest();
    int populationSize() const;
    Population evolvePopulation();
    Tour crossover(const Tour& parent1, const Tour& parent2);
    void mutate(Tour& tour);
    Tour tournamentSelection(Population& pop);
};

Population::Population(int populationSize, bool initialise) : tours(populationSize)
{
    if (initialise)
    {
        for (int i = 0; i < this->populationSize(); i++)
        {
            Tour newTour;
            newTour.generateIndividual();
            saveTour(i, newTour);
        }
    }
}
void Population::saveTour(int index, const Tour& tour)
{
    tours[index] = tour;
}


Tour Population::getTour(int index) const
{
    return tours[index];
}


Tour Population::getFittest()
{
    Tour fittest = tours[0];
    for (int i = 1; i < populationSize(); i++)
    {
        if (fittest.getFitness() <= getTour(i).getFitness())
        {
            fittest = getTour(i);
        }
    }
    return fittest;
}


int Population::populationSize() const
{
    return tours.capacity();
}

Population Population::evolvePopulation()
{
    Population newPopulation(this->populationSize(), false);

    for (int i = 0; i < newPopulation.populationSize(); i++)
    {
        Tour parent1 = tournamentSelection(*this);
        Tour parent2 = tournamentSelection(*this);
        Tour child = crossover(parent1, parent2);
        newPopulation.saveTour(i, child);
    }
    for (int i = 0; i < newPopulation.populationSize(); i++)
    {
        mutate((Tour&)(newPopulation.getTour(i)));
    }
    return newPopulation;
}

Tour Population::crossover(const Tour& parent1, const Tour& parent2)
{
    Tour child;
    int startPos = generateRandomInteger(parent1.tourSize());
    int endPos = generateRandomInteger(parent1.tourSize());
    for (int i = 0; i < child.tourSize(); i++)
    {
        if (startPos < endPos && i > startPos && i < endPos)
        {
            child.setCity(i, parent1.getCity(i));
        }
        else if (startPos > endPos)
        {
            if (!(i < startPos && i > endPos))
            {
                child.setCity(i, parent1.getCity(i));
            }
        }
    }
    for (int i = 0; i < parent2.tourSize(); i++)
    {
        if (!child.containsCity(parent2.getCity(i)))
        {
            // Loop to find a spare position in the child's tour
            for (int ii = 0; ii < child.tourSize(); ii++)
            {
                // Spare position found, add city
                if (child.getCity(ii).getX() == 0.0 && child.getCity(ii).getY() == 0.0)
                {
                    child.setCity(ii, parent2.getCity(i));
                    break;
                }
            }
        }
    }
    return child;
}

void Population::mutate(Tour& tour)
{
    for (int tourPos1 = 0; tourPos1 < tour.tourSize(); tourPos1++)
    {
        if (generateRandomDouble() < mutationRate)
        {
            int tourPos2 = generateRandomInteger(tour.tourSize());
            City city1 = tour.getCity(tourPos1);
            City city2 = tour.getCity(tourPos2);
            tour.setCity(tourPos2, city1);
            tour.setCity(tourPos1, city2);
        }
    }
}

Tour Population::tournamentSelection(Population& pop)
{
    Population tournament(tournamentSize, false);
    for (int i = 0; i < tournamentSize; i++)
    {
        int randomId = generateRandomInteger(pop.populationSize());
        tournament.saveTour(i, pop.getTour(randomId));
    }
    Tour fittest = tournament.getFittest();
    return fittest;
}


int main()
{
    destinationCities.push_back(City(0.000190032, -0.000285946));//Aberystwyth
    destinationCities.push_back(City(383.458, -0.000608756));//Brighton
    destinationCities.push_back(City(-27.0206, -282.758));//Edinburgh
    destinationCities.push_back(City(335.751, -269.577));//Exeter
    destinationCities.push_back(City(69.4331, -246.780));//Glasgow
    destinationCities.push_back(City(168.521, 31.4012));//Inverness
    destinationCities.push_back(City(320.350, -160.900));//Liverpool
    destinationCities.push_back(City(179.933, -318.031));//London
    destinationCities.push_back(City(492.671, -131.563));//Newcastle
    destinationCities.push_back(City(112.198, -110.561));//Nottingham
    destinationCities.push_back(City(306.320, -108.090));//Oxford
    destinationCities.push_back(City(217.343, -447.089));//Stratford
/*
    for (int i = 0; i < 20; i++)
    {
        destinationCities.push_back(City());
    }  */
    /*
        destinationCities.push_back(City(0.190032E-03, -0.285946E-03));//Aberystwyth
        destinationCities.push_back(City(168.521, 31.4012));//Inverness
        destinationCities.push_back(City(112.198, -110.561));//Nottingham
        destinationCities.push_back(City(69.4331, -246.780));//Glasgow
        destinationCities.push_back(City(-27.0206, -282.758));//Edinburgh
        destinationCities.push_back(City(179.933, -318.031));//London
        destinationCities.push_back(City(217.343, -447.089));//Stratford
        destinationCities.push_back(City(335.751, -269.577));//Exeter
        destinationCities.push_back(City(320.350, -160.900));//Liverpool
        destinationCities.push_back(City(306.320, -108.090));//Oxford
        destinationCities.push_back(City(383.458, -0.608756E-03));//Brighton
        destinationCities.push_back(City(492.671, -131.563));//Newcastle */


    Population pop(900, true);
    int initialDistance = pop.getFittest().getDistance();

    for (int i = 0; i < 35; i++)
    {
        cout << i << "th path: " << pop.getFittest().getDistance() << endl;
        pop = pop.evolvePopulation();
    }
    cout << "Solution:" << pop.getFittest() << endl;
    cout << "\nInfo:" << endl;
    cout << "Initial distance: " << initialDistance << endl;
    cout << "Solution distance: " << pop.getFittest().getDistance() << endl;
    return 0;
}
/*

Aberystwyth  0.190032E-03, -0.285946E-03
Brighton     383.458, -0.608756E-03
Edinburgh    - 27.0206, -282.758
Exeter       335.751, -269.577
Glasgow      69.4331, -246.780
Inverness    168.521, 31.4012
Liverpool    320.350, -160.900
London       179.933, -318.031
Newcastle    492.671, -131.563
Nottingham   112.198, -110.561
Oxford       306.320, -108.090
Stratford    217.343, -447.089


The optimal total distance is:
1595.738522033024
*/