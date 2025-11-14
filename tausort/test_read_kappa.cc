// Simple standalone test to verify reading kappa_5_band.dat
// Matches the layout used in RTS::load_bins up to reading tab_T and tab_p.

#include <iostream>
#include <fstream>
#include <vector>

int main(int argc, char **argv)
{
    const char *kap_name = (argc > 1) ? argv[1] : "kappa_5_band.dat";

    std::ifstream fp_rt(kap_name, std::ios::in | std::ios::binary);
    if (!fp_rt.is_open())
    {
        std::cerr << "Could not open file: " << kap_name << std::endl;
        return 1;
    }

    int N5000, NT, Np, Nbands;
    int pt_rhot, fullodf, scatter, back_heating;

    fp_rt.read(reinterpret_cast<char *>(&N5000), sizeof(int));
    fp_rt.read(reinterpret_cast<char *>(&NT), sizeof(int));
    fp_rt.read(reinterpret_cast<char *>(&Np), sizeof(int));
    fp_rt.read(reinterpret_cast<char *>(&Nbands), sizeof(int));
    fp_rt.read(reinterpret_cast<char *>(&pt_rhot), sizeof(int));
    fp_rt.read(reinterpret_cast<char *>(&fullodf), sizeof(int));
    fp_rt.read(reinterpret_cast<char *>(&scatter), sizeof(int));
    fp_rt.read(reinterpret_cast<char *>(&back_heating), sizeof(int));

    if (!fp_rt.good())
    {
        std::cerr << "Error while reading header from " << kap_name << std::endl;
        return 1;
    }

    std::cout << "Header values from " << kap_name << ":\n";
    std::cout << "  N5000          = " << N5000 << "\n";
    std::cout << "  NT             = " << NT << "\n";
    std::cout << "  Np             = " << Np << "\n";
    std::cout << "  Nbands         = " << Nbands << "\n";
    std::cout << "  pt_rhot        = " << pt_rhot << "\n";
    std::cout << "  fullodf        = " << fullodf << "\n";
    std::cout << "  scatter        = " << scatter << "\n";
    std::cout << "  back_heating   = " << back_heating << "\n";

    std::vector<double> tab_T(NT);
    std::vector<double> tab_p(Np);

    fp_rt.read(reinterpret_cast<char *>(tab_T.data()), NT * sizeof(double));
    fp_rt.read(reinterpret_cast<char *>(tab_p.data()), Np * sizeof(double));

    if (!fp_rt.good())
    {
        std::cerr << "Error while reading tab_T/tab_p from " << kap_name << std::endl;
        return 1;
    }

    std::cout.setf(std::ios::scientific);
    std::cout.precision(8);

    std::cout << "\nFirst few tab_T values (log10 T):\n";
    for (int i = 0; i < std::min(NT, 5); ++i)
        std::cout << "  tab_T[" << i << "] = " << tab_T[i] << "\n";

    std::cout << "\nFirst few tab_p values (log10 p):\n";
    for (int i = 0; i < std::min(Np, 5); ++i)
        std::cout << "  tab_p[" << i << "] = " << tab_p[i] << "\n";

    std::cout << "\nRead test finished successfully.\n";

    return 0;
}
