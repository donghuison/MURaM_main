/*
 * Simplified serial version of parts of RTS (radiative transfer solver).
 *
 * This file is intended purely for learning and experimentation:
 *  - CPU-only, single-process (no MPI Cartesian topology or communicators)
 *  - No ACCH memory manager or OpenACC pragmas
 *  - Focuses on:
 *      * Reading the kappa-bin binary file (same layout as RTS::load_bins)
 *      * Holding the main opacity / source-function tables on the CPU
 *      * Performing the core kappa/B interpolation for a single band
 *
 * What is deliberately omitted compared to RTS in rt.cc:
 *  - MPI topology, column communicators, and neighbour rank bookkeeping
 *  - GPU data regions and device/host transfers (ACCH::Malloc/Free/UpdateGPU)
 *  - Full multi-angle integration driver and flux / heating-rate computation
 *
 * How to read this file:
 *  - The structure of load_bins_serial mirrors RTS::load_bins in rt.cc
 *    (around lines 604–753), but uses std::vector and plain C++ I/O.
 *  - compute_band_opacities mirrors the OpenACC loop in rt.cc
 *    (around lines 956–985), but as simple nested for-loops.
 *
 * This code is not wired into the MURaM build; it can be compiled as a
 * standalone helper / example if desired.
 */

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{
    // These mirror the definitions in rt.h but are duplicated here to keep
    // this file independent of MPI / ACCH / OpenACC.
    static const int Nlam_serial = 328; // Number of intervals in ODF
    static const int Nbin_serial = 12;  // # bins per interval
    static const double TENLOG_serial = 2.30258509299405e0;
} // namespace

class RTS_Serial
{
public:
    // Publicly expose basic grid size and some fields for easy experimentation.
    int nx, ny, nz;

    std::vector<double> lgTe;   // log(T) per cell
    std::vector<double> lgPe;   // log(P) per cell
    std::vector<int> T_ind;     // temperature bin index per cell
    std::vector<int> P_ind;     // pressure bin index per cell
    std::vector<int> tr_switch; // transition-region switch (0/1) per cell

    std::vector<double> B;   // Planck function per cell (for a given band)
    std::vector<double> kap; // opacity per cell (for a given band)

    // Constructor: kap_name is the kappa-bin file; rttype/eps_const are used
    // only to reproduce some of the branching logic from RTS::load_bins.
    RTS_Serial(const std::string &kap_name,
               int nx_in, int ny_in, int nz_in,
               int rttype_in = 0,
               double eps_const_in = 0.0)
        : nx(nx_in), ny(ny_in), nz(nz_in),
          Nbands(0), NT(0), Np(0), N5000(0), fullodf(0),
          scatter(0), Npp(0), rttype(rttype_in), eps_const(eps_const_in)
    {
        const std::size_t ncell = static_cast<std::size_t>(nx) *
                                  static_cast<std::size_t>(ny) *
                                  static_cast<std::size_t>(nz);

        lgTe.assign(ncell, 0.0);
        lgPe.assign(ncell, 0.0);
        T_ind.assign(ncell, 0);
        P_ind.assign(ncell, 0);
        tr_switch.assign(ncell, 1);

        B.assign(ncell, 0.0);
        kap.assign(ncell, 0.0);

        load_bins_serial(kap_name);
    }

    // Perform kappa/B interpolation for a single band (CPU-only, serial).
    // This mirrors the OpenACC loop around lines 956–985 in rt.cc:
    //
    //  B[ind] = exp(xt*B_tab[band*NT+l+1]+(1.-xt)*B_tab[band*NT+l]);
    //  kap[ind] =
    //    exp(xt*(xp*kap_tab[(band*NT+l+1)*Np+m+1]+(1.-xp)*kap_tab[(band*NT+l+1)*Np+m])+
    //    (1.-xt)*(xp*kap_tab[(band*NT+l)*Np+m+1]+(1.-xp)*kap_tab[(band*NT+l)*Np+m]));
    //
    void compute_band_opacities(int band)
    {
        if (band < 0 || band >= Nbands)
        {
            std::cerr << "RTS_Serial::compute_band_opacities: band out of range\n";
            return;
        }
        if (NT < 2 || Np < 2)
        {
            std::cerr << "RTS_Serial::compute_band_opacities: NT/Np too small\n";
            return;
        }

        for (int y = 0; y < ny; ++y)
        {
            for (int x = 0; x < nx; ++x)
            {
                const int xyoff = y * nx * nz + x * nz;
                for (int z = 0; z < nz; ++z)
                {
                    const int ind = xyoff + z;

                    const int l = T_ind[ind];
                    const int m = P_ind[ind];

                    // These indices should satisfy 0 <= l <= NT-2, 0 <= m <= Np-2.
                    if (l < 0 || l >= NT - 1 || m < 0 || m >= Np - 1)
                    {
                        // For learning purposes, fail loudly if indices are out of range.
                        std::cerr << "RTS_Serial::compute_band_opacities: "
                                     "T_ind/P_ind out of range at cell "
                                  << ind << " (l=" << l << ", m=" << m << ")\n";
                        continue;
                    }

                    const double xt = (lgTe[ind] - tab_T[l]) * invT_tab[l];
                    const double xp = (lgPe[ind] - tab_p[m]) * invP_tab[m];

                    // Interpolate B_tab linearly in log-space.
                    const int B_offset = band * NT;
                    const double B_log =
                        xt * B_tab[B_offset + l + 1] +
                        (1.0 - xt) * B_tab[B_offset + l];
                    double B_val = std::exp(B_log);

                    // Interpolate kap_tab linearly in log-space.
                    const int base_up = (band * NT + l + 1) * Np;
                    const int base_down = (band * NT + l) * Np;

                    const double kap_up =
                        xp * kap_tab[base_up + m + 1] +
                        (1.0 - xp) * kap_tab[base_up + m];

                    const double kap_down =
                        xp * kap_tab[base_down + m + 1] +
                        (1.0 - xp) * kap_tab[base_down + m];

                    const double kap_log =
                        xt * kap_up +
                        (1.0 - xt) * kap_down;

                    double kap_val = std::exp(kap_log);

                    // TR-switch: zero out above the transition region.
                    const int sw = tr_switch[ind];
                    kap_val *= sw;
                    B_val *= sw;

                    kap[ind] = kap_val;
                    B[ind] = B_val;
                }
            }
        }
    }

private:
    // Kappa-table parameters (mirror rt.h for the subset we need).
    int Nbands, NT, Np;
    int N5000;
    int fullodf;
    int scatter;
    int Npp;

    int rttype;
    double eps_const;

    // band-P-T tabulated opacities and helpers.
    std::vector<double> tab_T;
    std::vector<double> invT_tab;
    std::vector<double> tab_p;
    std::vector<double> invP_tab;
    std::vector<float> kap_tab;

    // band-T tabulated source function.
    std::vector<float> B_tab;

    // Reference wavelength grid (optional).
    std::vector<float> B_5000_tab;
    std::vector<float> kap_5000_tab;

    // For full ODF (optional; stored as flat arrays).
    std::vector<float> nu_tab;
    std::vector<float> acont_pT; // size Nlam_serial*NT*Np when used
    std::vector<float> kcont_pT; // size Nlam_serial*NT*Np when used

    // Scattering tables for PP approximation (optional).
    std::vector<double> tau_pp_tab;
    std::vector<double> invtau_pp_tab;
    std::vector<float> kap_pp_tab; // size Nbands*Npp when used
    std::vector<float> abn_pp_tab; // size Nbands*Npp when used
    std::vector<float> sig_pp_tab; // size Nbands*Npp when used

    void load_bins_serial(const std::string &kap_name)
    {
        std::ifstream fp_rt(kap_name.c_str(),
                            std::ios::in | std::ios::binary);

        int pt_rhot = 0;
        int junk4 = 0;

        if (!fp_rt.is_open())
        {
            std::cerr << "RTS_Serial::load_bins_serial: kappa file not found: "
                      << kap_name << "\n";
            std::exit(1);
        }

        fp_rt.read(reinterpret_cast<char *>(&N5000), sizeof(int));
        fp_rt.read(reinterpret_cast<char *>(&NT), sizeof(int));
        fp_rt.read(reinterpret_cast<char *>(&Np), sizeof(int));
        fp_rt.read(reinterpret_cast<char *>(&Nbands), sizeof(int));
        fp_rt.read(reinterpret_cast<char *>(&pt_rhot), sizeof(int));
        fp_rt.read(reinterpret_cast<char *>(&fullodf), sizeof(int));
        fp_rt.read(reinterpret_cast<char *>(&scatter), sizeof(int));
        fp_rt.read(reinterpret_cast<char *>(&junk4), sizeof(int));

        if (pt_rhot != 0)
        {
            std::cerr << "RTS_Serial::load_bins_serial: pt_rhot = "
                      << pt_rhot
                      << ", but rho-T bins are not implemented in this serial version.\n";
            std::exit(1);
        }

        // This mirrors the original safety check: for scattering RT we need
        // either scattering bins or a photon-destruction probability.
        if ((scatter == 0) && (eps_const == 0.0) && (rttype == 1))
        {
            std::cerr << "RTS_Serial::load_bins_serial: rt_type = "
                      << rttype << ", scatter = " << scatter
                      << ", eps_const = " << eps_const
                      << " -> invalid combination for scattering.\n";
            std::exit(1);
        }

        tab_T.assign(NT, 0.0);
        tab_p.assign(Np, 0.0);
        invT_tab.assign(NT, 0.0);
        invP_tab.assign(Np, 0.0);

        fp_rt.read(reinterpret_cast<char *>(tab_T.data()),
                   static_cast<std::streamsize>(NT * sizeof(double)));
        fp_rt.read(reinterpret_cast<char *>(tab_p.data()),
                   static_cast<std::streamsize>(Np * sizeof(double)));

        for (int i = 0; i < NT; ++i)
            tab_T[i] *= TENLOG_serial;
        for (int j = 0; j < Np; ++j)
            tab_p[j] *= TENLOG_serial;

        for (int l = 0; l <= NT - 2; ++l)
            invT_tab[l] = 1.0 / (tab_T[l + 1] - tab_T[l]);
        for (int m = 0; m <= Np - 2; ++m)
            invP_tab[m] = 1.0 / (tab_p[m + 1] - tab_p[m]);

        if (N5000)
        {
            kap_5000_tab.assign(static_cast<std::size_t>(NT) *
                                    static_cast<std::size_t>(Np),
                                0.0f);
            B_5000_tab.assign(NT, 0.0f);

            fp_rt.read(reinterpret_cast<char *>(kap_5000_tab.data()),
                       static_cast<std::streamsize>(NT * Np * sizeof(float)));
            fp_rt.read(reinterpret_cast<char *>(B_5000_tab.data()),
                       static_cast<std::streamsize>(NT * sizeof(float)));
        }

        kap_tab.assign(static_cast<std::size_t>(Nbands) *
                           static_cast<std::size_t>(NT) *
                           static_cast<std::size_t>(Np),
                       0.0f);
        B_tab.assign(static_cast<std::size_t>(Nbands) *
                         static_cast<std::size_t>(NT),
                     0.0f);

        fp_rt.read(reinterpret_cast<char *>(kap_tab.data()),
                   static_cast<std::streamsize>(kap_tab.size() * sizeof(float)));
        fp_rt.read(reinterpret_cast<char *>(B_tab.data()),
                   static_cast<std::streamsize>(B_tab.size() * sizeof(float)));

        if (fullodf)
        {
            // Full ODF mode: read nu_tab and continuum opacity tables.
            nu_tab.assign(Nlam_serial, 0.0f);
            fp_rt.read(reinterpret_cast<char *>(nu_tab.data()),
                       static_cast<std::streamsize>(Nlam_serial * sizeof(float)));

            if (scatter > 0)
            {
                const std::size_t odf_size =
                    static_cast<std::size_t>(Nlam_serial) *
                    static_cast<std::size_t>(NT) *
                    static_cast<std::size_t>(Np);

                acont_pT.assign(odf_size, 0.0f);
                kcont_pT.assign(odf_size, 0.0f);

                fp_rt.read(reinterpret_cast<char *>(acont_pT.data()),
                           static_cast<std::streamsize>(odf_size * sizeof(float)));
                fp_rt.read(reinterpret_cast<char *>(kcont_pT.data()),
                           static_cast<std::streamsize>(odf_size * sizeof(float)));

                Npp = 1;
                tau_pp_tab.assign(Npp, 0.0);
                invtau_pp_tab.assign(Npp, 0.0);

                tau_pp_tab[0] = -99.0;
                invtau_pp_tab[0] = 1.0;

                // If we have scattering bins and full ODF but run non-scattering RT
                // (rttype == 0), we need to add acont to kappa in log-space.
                if (rttype == 0)
                {
                    for (int lam = 0; lam < Nlam_serial; ++lam)
                    {
                        for (int bin = 0; bin < Nbin_serial; ++bin)
                        {
                            for (int tt = 0; tt < NT; ++tt)
                            {
                                for (int pp = 0; pp < Np; ++pp)
                                {
                                    const int band_idx = Nbin_serial * lam + bin;
                                    const std::size_t kap_index =
                                        (static_cast<std::size_t>(band_idx) * NT + tt) * Np + pp;
                                    const std::size_t odf_index =
                                        (static_cast<std::size_t>(lam) * NT + tt) * Np + pp;

                                    const double kappa_val =
                                        std::exp(static_cast<double>(kap_tab[kap_index]));
                                    const double acont_val =
                                        static_cast<double>(acont_pT[odf_index]);

                                    kap_tab[kap_index] =
                                        static_cast<float>(std::log(kappa_val + acont_val));
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (scatter > 0)
        {
            // Non-ODF scattering case: read band-dependent scattering tables.
            Npp = scatter;
            scatter = 1;

            const std::size_t sig_abn_size =
                static_cast<std::size_t>(Nbands) *
                static_cast<std::size_t>(NT) *
                static_cast<std::size_t>(Np);

            abn_tab_flat.assign(sig_abn_size, 0.0f);
            sig_tab_flat.assign(sig_abn_size, 0.0f);

            fp_rt.read(reinterpret_cast<char *>(abn_tab_flat.data()),
                       static_cast<std::streamsize>(sig_abn_size * sizeof(float)));
            fp_rt.read(reinterpret_cast<char *>(sig_tab_flat.data()),
                       static_cast<std::streamsize>(sig_abn_size * sizeof(float)));

            tau_pp_tab.assign(Npp, 0.0);
            invtau_pp_tab.assign(Npp, 0.0);

            fp_rt.read(reinterpret_cast<char *>(tau_pp_tab.data()),
                       static_cast<std::streamsize>(Npp * sizeof(double)));

            for (int z = 0; z <= Npp - 2; ++z)
                invtau_pp_tab[z] =
                    1.0 / (tau_pp_tab[z + 1] - tau_pp_tab[z]);

            const std::size_t pp_size =
                static_cast<std::size_t>(Nbands) *
                static_cast<std::size_t>(Npp);

            kap_pp_tab.assign(pp_size, 0.0f);
            abn_pp_tab.assign(pp_size, 0.0f);
            sig_pp_tab.assign(pp_size, 0.0f);

            fp_rt.read(reinterpret_cast<char *>(abn_pp_tab.data()),
                       static_cast<std::streamsize>(pp_size * sizeof(float)));
            fp_rt.read(reinterpret_cast<char *>(kap_pp_tab.data()),
                       static_cast<std::streamsize>(pp_size * sizeof(float)));
            fp_rt.read(reinterpret_cast<char *>(sig_pp_tab.data()),
                       static_cast<std::streamsize>(pp_size * sizeof(float)));
        }

        if ((scatter == 0) && (rttype == 1))
        {
            // Scattering RT but no scatter bins: create a trivial tau-grid.
            Npp = 1;
            tau_pp_tab.assign(Npp, 0.0);
            invtau_pp_tab.assign(Npp, 0.0);

            tau_pp_tab[0] = -99.0;
            invtau_pp_tab[0] = 1.0;
        }

        fp_rt.close();
    }

    // Flat storage for abn_tab and sig_tab from the scattering branch.
    std::vector<float> abn_tab_flat;
    std::vector<float> sig_tab_flat;
};

#ifdef RT_SERIAL_EXAMPLE
// Simple example / smoke test (not integrated into the main code):
//  - Adjust the file name below to point to a valid kappa-bin file.
int main(int argc, char **argv)
{
    const std::string kap_name =
        (argc > 1 ? std::string(argv[1]) : std::string("kappa_bins.dat"));

    const int nx = 4, ny = 4, nz = 4;

    RTS_Serial rts_serial(kap_name, nx, ny, nz, /*rttype=*/0, /*eps_const=*/0.0);

    const std::size_t ncell =
        static_cast<std::size_t>(nx) *
        static_cast<std::size_t>(ny) *
        static_cast<std::size_t>(nz);

    // Fill some trivial test data for lgTe, lgPe, indices, and tr_switch.
    for (std::size_t i = 0; i < ncell; ++i)
    {
        rts_serial.lgTe[i] = std::log(5777.0);
        rts_serial.lgPe[i] = std::log(1.0e5);
        rts_serial.T_ind[i] = 0;
        rts_serial.P_ind[i] = 0;
        rts_serial.tr_switch[i] = 1;
    }

    // Compute opacities for band 0 and print a few sample values.
    rts_serial.compute_band_opacities(0);

    std::cout << "Sample B and kap values (first few cells):\n";
    for (int i = 0; i < 5 && i < static_cast<int>(ncell); ++i)
    {
        std::cout << "cell " << i
                  << " B=" << rts_serial.B[i]
                  << " kap=" << rts_serial.kap[i] << "\n";
    }

    return 0;
}
#endif
