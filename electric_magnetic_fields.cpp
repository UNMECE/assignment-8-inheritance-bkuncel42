#include <iostream>
#include <cmath>

const double EPSILON_0 = 8.854e-12;
const double MU_0 = 4 * M_PI * 1e-7;

class Field {
protected:
    double* components; // x, y, z

public:
    Field() {
        components = new double[3]{0.0, 0.0, 0.0};
    }

    Field(double x, double y, double z) {
        components = new double[3]{x, y, z};
    }

    Field(const Field& other) {
        components = new double[3];
        for (int i = 0; i < 3; ++i)
            components[i] = other.components[i];
    }

    Field& operator=(const Field& other) {
        if (this != &other) {
            for (int i = 0; i < 3; ++i)
                components[i] = other.components[i];
        }
        return *this;
    }

    virtual ~Field() {
        delete[] components;
    }

    virtual void printMagnitude() const {
        std::cout << "Field components: (" << components[0] << ", " << components[1] << ", " << components[2] << ")\n";
    }

    double magnitude() const {
        return std::sqrt(components[0]*components[0] + components[1]*components[1] + components[2]*components[2]);
    }

    double getX() const { return components[0]; }
    double getY() const { return components[1]; }
    double getZ() const { return components[2]; }
};

// ---------------- Electric Field ----------------
class Electric_Field : public Field {
private:
    double calculated_E;

public:
    Electric_Field() : Field(), calculated_E(0) {}
    Electric_Field(double x, double y, double z) : Field(x, y, z), calculated_E(0) {}

    Electric_Field(const Electric_Field& other) : Field(other), calculated_E(other.calculated_E) {}

    void calculateElectricField(double Q, double r) {
        calculated_E = Q / (4 * M_PI * EPSILON_0 * r * r);
    }

    void printCalculatedElectricField() const {
        std::cout << "Calculated Electric Field Magnitude: " << calculated_E << " N/C\n";
    }

    Electric_Field operator+(const Electric_Field& other) const {
        return Electric_Field(
            components[0] + other.components[0],
            components[1] + other.components[1],
            components[2] + other.components[2]
        );
    }

    friend std::ostream& operator<<(std::ostream& os, const Electric_Field& ef) {
        os << "Electric Field: (" << ef.components[0] << ", " << ef.components[1] << ", " << ef.components[2] << ")";
        return os;
    }
};

// ---------------- Magnetic Field ----------------
class Magnetic_Field : public Field {
private:
    double calculated_B;

public:
    Magnetic_Field() : Field(), calculated_B(0) {}
    Magnetic_Field(double x, double y, double z) : Field(x, y, z), calculated_B(0) {}

    Magnetic_Field(const Magnetic_Field& other) : Field(other), calculated_B(other.calculated_B) {}

    void calculateMagneticField(double I, double r) {
        calculated_B = MU_0 * I / (2 * M_PI * r);
    }

    void printCalculatedMagneticField() const {
        std::cout << "Calculated Magnetic Field Magnitude: " << calculated_B << " T\n";
    }

    Magnetic_Field operator+(const Magnetic_Field& other) const {
        return Magnetic_Field(
            components[0] + other.components[0],
            components[1] + other.components[1],
            components[2] + other.components[2]
        );
    }

    void calculateUnitVector(double& ux, double& uy, double& uz) const {
        double mag = magnitude();
        if (mag == 0) {
            ux = uy = uz = 0.0;
        } else {
            ux = components[0] / mag;
            uy = components[1] / mag;
            uz = components[2] / mag;
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const Magnetic_Field& mf) {
        os << "Magnetic Field: (" << mf.components[0] << ", " << mf.components[1] << ", " << mf.components[2] << ")";
        return os;
    }
};

// ---------------- Main ----------------
int main() {
    // Create Electric Fields
    Electric_Field e1(1e5, 10.9, 1.7e2);
    Electric_Field e2(300, 400, 1200);

    std::cout << e1 << "\n";
    e1.printMagnitude();
    std::cout << "Magnitude: " << e1.magnitude() << "\n";

    std::cout << e2 << "\n";
    std::cout << "Magnitude: " << e2.magnitude() << "\n";

    e1.calculateElectricField(1e-6, 0.05);
    e1.printCalculatedElectricField();

    // Create Magnetic Fields
    Magnetic_Field m1(0.5, 1.5, 2.5);
    Magnetic_Field m2(1.0, 1.0, 1.0);

    std::cout << m1 << "\n";
    m1.printMagnitude();
    std::cout << "Magnitude: " << m1.magnitude() << "\n";

    m1.calculateMagneticField(5.0, 0.01);
    m1.printCalculatedMagneticField();

    double ux, uy, uz;
    m1.calculateUnitVector(ux, uy, uz);
    std::cout << "Magnetic Unit Vector: (" << ux << ", " << uy << ", " << uz << ")\n";

    // Operator Overload Demo
    Electric_Field e3 = e1 + e2;
    Magnetic_Field m3 = m1 + m2;

    std::cout << "\nResult of e1 + e2: " << e3 << "\n";
    std::cout << "Result of m1 + m2: " << m3 << "\n";

    return 0;
}

