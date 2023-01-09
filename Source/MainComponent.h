#pragma once

#include <JuceHeader.h>

//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent  : public juce::AudioAppComponent,
                       public juce::Slider::Listener
{
public:
    //==============================================================================
    MainComponent();
    ~MainComponent() override;

    void mouseDrag(const MouseEvent& e) override;
    void mouseDown(const MouseEvent& e) override;
    void mouseUp(const MouseEvent& e) override;
    ////void mouseExit(const MouseEvent& e) override;


    //==============================================================================
    void prepareToPlay(int samplesPerBlockExpected, double sampleRate) override;
    void getNextAudioBlock(const juce::AudioSourceChannelInfo& bufferToFill) override;
    void releaseResources() override;

    //==============================================================================
    void paint(juce::Graphics& g) override;
    void resized() override;
    void sliderValueChanged(Slider* slider) override;


private:
    //==============================================================================
    // Your private member variables go here...

    Slider frequencySlider;
    Label frequencyLabel;

    Slider sig_0_Slider;
    Label sig_0_Label;

    Slider sig_1_Slider;
    Label sig_1_Label;

    // GENERAL PARAMETERS
    double fs;
    //float gain = 1000000000;
    float gain = 1000;
    float t = 0;
    float preOutput;
    double k; // time step in Finite Difference scheme


    // BOWING MODEL PARAMETERS
    double mu_C = 0.3;
    double mu_S = 0.8;
    double fN = 10;
    //double fC = fN * mu_C;
    //double fS = fN * mu_S;
    double fC;
    double fS;
    double vB = 0.1;
    double vS = 0.1;
    double s0 = 1e5;
    double s1 = 0.001 * sqrt(s0);
    //double s1 = 0.0;
    double s2 = 4;
    double s3_fac = 0;
    double s3;
    double w_rnd;
    double z_ba;
    double tol = 0.0000001;

    double x_bow_membrane = 0.2; // absolute
    double y_bow_membrane = 0.125; // absolute

    double f0;

    // MEMBRANE PARAMETERS
    double L_membrane = 0.3; // working with real dimensions
    double rho_membrane = 1400; // Density of mylar - what some drum heads are made of..
    double H_membrane = 0.007; // [m] : check: https://www.andertons.co.uk/drum-skins-drum-heads-guide
    double T_membrane = 10000; // [N/m]

    //// Add these as T60 stuff.. 
    //double* lossFreqs_membrane = new double[2]{ 100, 3000 }; 
    //double* lossT60s_membrane = new double[2]{ 3, 1 };
    //double zeta1_plate;
    //double zeta2_plate;

    double sig0;
    double sig1;

    double c_membrane;

    double h_membrane; // spatial step in FD scheme
    double N_membrane; // no of spatial points to discretize string 
    int Nx_membrane;
    int Ny_membrane;

    // MEMBRANE CONNECTION:
    // there is no single connection point for the membrane.. just a distrib. 
    std::vector<double> I_M; // interpolant grid
    std::vector<double> J_M; // spreading function grid

    std::vector<std::vector<double>> locsDo;


    // TUBE PARAMETERS
    double rho_tube = 1.225; // [kg / m ^ 3]
    double radius = L_membrane * 0.5;
    double A_tube = double_Pi * radius * radius; // [m ^ 2]
    double T_tube = 80; // [N] [kg*m/s^2]
    double L_tube = 0.4;

    double c_tube;

    double h_tube; // spatial step in FD scheme
    double N_tube; // no of spatial points to discretize string 

    double sig0_tube = 0; // no damping along tube yet.
    double sig1_tube = 0; // no damping along tube yet.

    // TUBE RADIATION END PARAMS:
    double alpha1;
    double alpha2;

    double F1;
    double F2;
    double F3;

    // TUBE CONNECTION:
    // there is no single connection point for the tube.. just a distrib. 
    std::vector<double> I_T; // interpolant grid
    std::vector<double> J_T; // spreading function grid


    //double freq = 300; // initial string freq
    //double c = sqrt(T / (rho * H));
    //double gamma = freq * 2 * L; // gamma param in the PDE (wave speed)
    //double lossFreqs = [100, 4; 5000, 2]; // w1, T60(w1) w2, T60(w2)


    //double* lossFreqs = new double[2]{ 100, 5000 }; 
    //double* lossT60s = new double[2]{ 4, 2 }; 

    // PLATE OUTPUT:
    double x_plate_out = 0.43; // absolute value
    double y_plate_out = 0.15; // absolute value


    // BOWING POSITION INIT

    /// Parameters you may want to modulate:
    double bp = 0.29; // in percentage
    int bP;
    double alpha_bow;

    std::vector<double> I_B; // interpolant grid for bowing pos
    std::vector<double> J_B; // spreading function grid for bowing pos

    int l_inp;
    int l_inp_plus1;
    int m_inp;
    int m_inp_plus1;
    double alpha_x_inp;
    double alpha_y_inp;



    // Constants for update eq
    double C1;
    double C2;
    double C3;
    double C4;

    /// end global parameters




    // pointers to MEMBRANE states
    std::vector<double*> u;

    // states
    std::vector<std::vector<double>> uVecs;

    double* uTmp = nullptr;


    // pointers to TUBE states
    std::vector<double*> w;

    // states
    std::vector<std::vector<double>> wVecs;

    double* wTmp = nullptr;


    // INITIALIZE PARAMS FOR NEWTON-RAPHSON

    double z_prev;
    double r_last;
    double r_prev;
    double a_last;
    double a_prev;

    double vrel_last;
    double z_last;

    double F_last;

    // NEWTON-RAPHSON SOLVER
    double z_ss;
    double alpha_fct;
    double g1;
    double g2;
    double z_ss_d_vrel;
    double alpha_fct_d_vrel;
    double alpha_fct_d_z;
    double r_last_d_vrel;
    double r_last_d_z;
    double g1_d_vrel;
    double g1_d_z;
    double g2_d_vrel;
    double g2_d_z;
    double fac_1_over_det;
    double A_mat_new;
    double B_mat_new;
    double C_mat_new;
    double D_mat_new;
    double vrel;
    double z;
    double r;
    double a;
    double Fbow;

    // CONNECTION FORCE
    double F;


    //std::vector<double> uPrev; // prev values of string displacement
    //std::vector<double> u; // prev values of string displacement
    //std::vector<double> uNext; // prev values of string displacement
    //std::vector<double> I_grid; // interpolant grid
    //std::vector<double> J_grid; // spreading function grid

    //std::vector<double> output_interim; 

    //std::vector<double*> I_grid;
    //std::vector<double*> J_grid; // spreading function grid


    /// buffer related parameters
    std::vector<double> myBuffer; // am I using this.. no ?

    // Mouse control parameters
    float opacity = 0.0;
    int mouseX;
    int mouseY;
    int mouseX_up;
    int mouseY_up;
    int mouseX_down;
    int mouseY_down;
    //float FB_var;
    double fN_var = fN;
    double vB_var = vB;
    double fC_var;
    double fS_var;
    double z_ba_var;
    double s3_var;
    float x_inp_var;
    float y_inp_var;
    float x_out_var;
    float y_out_var;
    float x_var = 0;
    float y_var = 0;
    float xPrev = 0;
    float yPrev = 0;
    float resultant = 0;
    float resultantPrev = 0;
    float vel;


    std::deque<double> x_positions;
    std::deque<double> y_positions;



    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(MainComponent)
};
