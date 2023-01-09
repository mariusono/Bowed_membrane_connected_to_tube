#include "MainComponent.h"

//==============================================================================
MainComponent::MainComponent()
{
    // Make sure you set the size of the component after
    // you add any child components.
    setSize(800, 600);

    bool val = (juce::RuntimePermissions::isRequired(juce::RuntimePermissions::recordAudio)
        && !juce::RuntimePermissions::isGranted(juce::RuntimePermissions::recordAudio));

    // Some platforms require permissions to open input channels so request that here
    if (juce::RuntimePermissions::isRequired(juce::RuntimePermissions::recordAudio)
        && !juce::RuntimePermissions::isGranted(juce::RuntimePermissions::recordAudio))
    {
        juce::RuntimePermissions::request(juce::RuntimePermissions::recordAudio,
            [&](bool granted) { setAudioChannels(granted ? 2 : 0, 2); });
    }
    else
    {

        // Specify the number of input and output channels that we want to open
        setAudioChannels(2, 2);
    }



    // Sliders !
    addAndMakeVisible(frequencySlider);
    frequencySlider.setRange(30.0, 650.0);          // [1]
    frequencySlider.setTextValueSuffix(" Hz");     // [2]
    frequencySlider.addListener(this);             // [3]

    addAndMakeVisible(frequencyLabel);
    frequencyLabel.setText("Frequency", juce::dontSendNotification);
    frequencyLabel.attachToComponent(&frequencySlider, true); // [4]

    frequencySlider.setValue(220.0); // [5]

    addAndMakeVisible(sig_0_Slider);
    sig_0_Slider.setRange(0.0, 10.0);          // [1]
    sig_0_Slider.setTextValueSuffix(" [-]");     // [2]
    sig_0_Slider.addListener(this);             // [3]

    addAndMakeVisible(sig_0_Label);
    sig_0_Label.setText("Freq Independent Damping", juce::dontSendNotification);
    sig_0_Label.attachToComponent(&sig_0_Slider, true); // [4]

    sig_0_Slider.setValue(0.0); // [5]

    addAndMakeVisible(sig_1_Slider);
    sig_1_Slider.setRange(0.0, 0.005);          // [1]
    sig_1_Slider.setTextValueSuffix(" [-]");     // [2]
    sig_1_Slider.addListener(this);             // [3]

    addAndMakeVisible(sig_1_Label);
    sig_1_Label.setText("Freq Independent Damping", juce::dontSendNotification);
    sig_1_Label.attachToComponent(&sig_1_Slider, true); // [4]

    sig_1_Slider.setValue(0.0); // [5]


}

MainComponent::~MainComponent()
{
    // This shuts down the audio device and clears the audio source.
    shutdownAudio();
}


//==============================================================================
// Other functions

double linearMapping(float rangeIn_top, float rangeIn_bottom, float rangeOut_top, float rangeOut_bottom, float value) {
    double newValue = rangeOut_bottom + ((rangeOut_top - rangeOut_bottom) * (value - rangeIn_bottom) / (rangeIn_top - rangeIn_bottom));
    return newValue;
}

template <typename T> int sign(T val) { // I'm not super sure what this template thing does.. just took it from a stackoverflow thread.
    return (T(0) < val) - (val < T(0));
}


std::vector<double> hamming(int N) { // Hamming distribution vector
    double alpha = 0.54;
    double beta = 1 - alpha;

    std::vector<double> h;
    h.resize(N);

    for (int i = 0; i < N; ++i) {
        h[i] = alpha - beta * cos(2 * double_Pi * i / (N - 1));
    }

    return h;
}

std::vector<double> hann(int N) { // Hann distribution vector
    double alpha = 0.5;
    double beta = 1 - alpha;

    std::vector<double> h;
    h.resize(N);

    for (int i = 0; i < N; ++i) {
        h[i] = alpha - beta * cos(2 * double_Pi * i / (N - 1));
    }

    return h;
}



//==============================================================================
// Mouse functions:
void MainComponent::mouseDrag(const MouseEvent& e)
{
    mouseX = e.x;
    mouseY = e.y;

    int timeVal = e.eventTime.getMilliseconds();
    //Logger::getCurrentLogger()->outputDebugString("timeVal: (" + juce::String(timeVal) + ")");

    opacity = 0.7;

    y_inp_var = (-(e.y / static_cast<float>(getHeight()) - 1));
    //y_inp_var = linearMapping(1, 0, 2.0, 0.0001, y_inp_var);

    fN_var = y_inp_var * fN;
    vB_var = y_inp_var * vB;

    fS_var = fN_var * mu_S;
    fC_var = fN_var * mu_C;
    z_ba_var = 0.7 * (mu_C * fN_var) / s0;
    s3_var = s3_fac * fN;


    x_inp_var = e.x / static_cast<float>(getWidth());
    x_inp_var = linearMapping(1, 0, 0.8, 0.2, x_inp_var);


 /*   resultant = sqrt((y_var - yPrev) * (y_var - yPrev) - (x_var - xPrev) * (x_var - xPrev));

    x_positions.pop_front();
    x_positions.push_back(x_inp_var);

    y_positions.pop_front();
    y_positions.push_back(y_inp_var);

    vel = (resultant - resultantPrev);*/

    repaint();


}

void MainComponent::mouseUp(const MouseEvent& e)
{
    mouseX_up = e.x;
    mouseY_up = e.y;

    opacity = 0;

    //Logger::getCurrentLogger()->outputDebugString("Mouse up at: (" + String(e.x) + ")");

    //// This stuff below can be put into a function called reset() .. somewhere

    //fN_var = fN;
    //vB_var = vB;
    //fN_var = 0.1; // if fN is 0, don't do the NR ! just add condition in the while of the NR
    fN_var = 0.0; // if fN is 0, don't do the NR ! just add condition in the while of the NR
    vB_var = 0.0;
    fS_var = fN_var * mu_S;
    fC_var = fN_var * mu_C;
    z_ba_var = 0.7 * (mu_C * fN_var) / s0;
    s3_var = s3_fac * fN;

    repaint();
}


void MainComponent::mouseDown(const MouseEvent& e) // not sure this is needed...
{

    fN_var = fN_var;
    vB_var = vB_var;
    fS_var = fN_var * mu_S;
    fC_var = fN_var * mu_C;
    z_ba_var = 0.7 * (mu_C * fN_var) / s0;
    s3_var = s3_fac * fN;

    //Logger::getCurrentLogger()->outputDebugString("Mouse down at: (" + String(e.x) + ")");

    repaint();
}



void MainComponent::sliderValueChanged(Slider* slider)
{

    if (slider == &frequencySlider)
    {
        f0 = frequencySlider.getValue();
        c_membrane = 2 * f0 * L_membrane;
    }
    else if (slider == &sig_0_Slider)
    {
        sig0 = sig_0_Slider.getValue();
        //Logger::getCurrentLogger()->outputDebugString(String(sig0));
    }
    else if (slider == &sig_1_Slider)
    {
        sig1 = sig_1_Slider.getValue();
        //Logger::getCurrentLogger()->outputDebugString(String(sig1));
    }

}





//==============================================================================
void MainComponent::prepareToPlay(int samplesPerBlockExpected, double sampleRate) // how do I change the sampleRate .. ?
{
    // This function will be called when the audio device is started, or when
    // its settings (i.e. sample rate, block size, etc) are changed.

    // You can use this function to initialise any resources you might need,
    // but be careful - it will be called on the audio thread, not the GUI thread.

    // For more details, see the help for AudioProcessor::prepareToPlay()

    // Sampling rate and time step
    fs = sampleRate;
    k = 1 / fs;


    // Print some info related to sample rate and audio block size
    Logger::getCurrentLogger()->outputDebugString("sampleRate at beginning of prepareToPlay: (" + String(sampleRate) + ")");
    Logger::getCurrentLogger()->outputDebugString("samplesPerBlockExpected at beginning of prepareToPlay: (" + String(samplesPerBlockExpected) + ")");
    Logger::getCurrentLogger()->outputDebugString("fs: (" + String(fs) + ")");


    // MEMBRANE STUFF: (this can be ported to a Membrane class..)
    
    // Can tweak these with the slider..
    sig0 = 4.5949250;
    sig1 = 0.0000264807;

    c_membrane = sqrt(T_membrane / (rho_membrane * H_membrane));

    // Stability condition. Revisit where this comes from. 
    h_membrane = sqrt(c_membrane * c_membrane * k * k + 4 * sig1 * k + sqrt(pow((c_membrane * c_membrane * k * k + 4 * sig1 * k), 2)));

    // Max no of elems
    N_membrane = floor(L_membrane / h_membrane);
    if (N_membrane > 18)
    {
        N_membrane = 18;
    }

    // Recalculate h_membrane:
    h_membrane = L_membrane / N_membrane;

    //  CIRCULAR GRID SELECTION 
    std::vector<double> xVec(N_membrane * N_membrane); // what is the difference to xVec(N_membrane,0); ?
    std::vector<double> yVec(N_membrane * N_membrane);
    std::vector<double> selLoc(2);
    double radius_circ = N_membrane / 2 - 2.5;
    double circ_center_x = floor(N_membrane * 0.5);
    double circ_center_y = floor(N_membrane * 0.5);
    //locsDo.resize(N_membrane * N_membrane,0);
    //v1.push_back(0);

    int idx_x_y;
    for (int iX = 0; iX < N_membrane; ++iX)
    {
        for (int iY = 0; iY < N_membrane; ++iY)
        {
            idx_x_y = (iY)+(iX)*N_membrane;
            xVec[idx_x_y] = iX;
            yVec[idx_x_y] = iY;
        }
    }


    for (int iX = 0; iX < N_membrane; ++iX)
    {
        for (int iY = 0; iY < N_membrane; ++iY)
        {
            idx_x_y = (iY)+(iX)*N_membrane;

            if ((xVec[idx_x_y] - circ_center_x) * (xVec[idx_x_y] - circ_center_x) + (yVec[idx_x_y] - circ_center_y) * (yVec[idx_x_y] - circ_center_y) < radius_circ * radius_circ) {
                selLoc[0] = iX;
                selLoc[1] = iY;

                locsDo.push_back(selLoc);
            }
        }
    }

    Logger::getCurrentLogger()->outputDebugString("locsDo size: (" + String(locsDo.size()) + ")");

    // TUBE STUFF:  (this can be ported to a Tube class..)

    c_tube = sqrt(T_tube / (rho_tube * A_tube));

    // Stability condition, no damping along the tube (only radiation at the tube end..):
    h_tube = c_tube * k;
    N_tube = floor(L_tube / h_tube);
    if (N_tube > 200)
    {
        N_tube = 200;
    }
    // Recalcualte h_tube:
    h_tube = L_tube / N_tube;


    // Radiation damping params:
    alpha1 = 1.0 / (4 * 0.6133 * 0.6133 * (c_tube / L_tube));
    alpha2 = L_tube / (0.6133 * sqrt(A_tube / double_Pi));

    // Some constants for tube update equation at radiation end:
    F1 = (k * k * c_tube * c_tube / (h_tube * h_tube));
    F2 = 2 * (1 - F1);
    F3 = F1;


    Logger::getCurrentLogger()->outputDebugString("N_membrane: (" + String(N_membrane) + ")");
    Logger::getCurrentLogger()->outputDebugString("N_tube: (" + String(N_tube) + ")");

    // DISTRIBUTIONS FOR CONNECTIONS

    // MEMBRANE:
    I_M.resize(N_membrane * N_membrane, 0);
    J_M.resize(N_membrane * N_membrane, 0);

    // Calculating connection distribution such that it is spread over half the area (centered) with a Hann distribution.
    std::vector<double> H1 = hann((int)floor(N_membrane/2));  // you could use hamming instead.. or try a uniform distrib.. 

    std::vector<double> H2;
    H2.resize(N_membrane);

    for (int i = floor(N_membrane / 4); i < floor(N_membrane / 2) + floor(N_membrane / 4); ++i) {
        H2[i] = H1[i- floor(N_membrane / 4)];
    }

    //int idx_x_y;
    for (int iX = 0; iX < N_membrane; ++iX)
    {
        for (int iY = 0; iY < N_membrane; ++iY)
        {
            idx_x_y = (iY) + (iX) * N_membrane;
            I_M[idx_x_y] = H2[iX] * H2[iY];
        }
    }

    // calculate sum of distribution (for normalization)
    double sum_I_M = 0;
    for (int j = 0; j < N_membrane * N_membrane; ++j) {
        sum_I_M = sum_I_M + I_M[j];
    }

    // normalize I_M such that the sum of all its elements gives 1 (so that it is a distribution !)
    for (int j = 0; j < N_membrane * N_membrane; ++j) {
        I_M[j] = I_M[j] / sum_I_M;
    }

    // Spreading function:
    for (int i = 0; i < I_M.size(); ++i)
    {
        J_M[i] = I_M[i] * (1 / (h_membrane * h_membrane)); // speed up: keep divisions out of loop ! 
    }


    // TUBE:
    I_T.resize(N_tube, 0);
    J_T.resize(N_tube, 0);

    std::vector<double> win_hann = hann((int) N_tube * 0.1);  // hann window over 10% of the length of the tube.. results in a connection over half this length

    for (int j = 0; j < floor(win_hann.size() / 2); ++j) {
        I_T[j] = win_hann[j + floor(win_hann.size() / 2)];
    }

    double sum_I_T = 0;
    for (int j = 0; j < N_tube; ++j) {
        sum_I_T = sum_I_T + I_T[j];
    }

    // normalize I_T such that the sum of all its elements gives 1 (so that it is a distribution !)
    for (int j = 0; j < N_tube; ++j) {
        I_T[j] = I_T[j] / sum_I_M;
    }

    // Spreading function
    for (int i = 0; i < I_T.size(); ++i)
    {
        J_T[i] = I_T[i] * (1 / h_tube); // speed up: keep divisions out of loop ! 
    }



    // Bow interpolation and spreading function init:
    I_B.resize(N_membrane * N_membrane, 0);
    J_B.resize(N_membrane * N_membrane, 0);


    // INITIALIZE MEMBRANE STATE VECTORS
    uVecs.reserve(3);

    for (int i = 0; i < 3; ++i)
        uVecs.push_back(std::vector<double>(N_membrane * N_membrane, 0));

    u.resize(3);

    for (int i = 0; i < u.size(); ++i)
        u[i] = &uVecs[i][0];



    // INITIALIZE TUBE STATE VECTORS
    wVecs.reserve(3);

    for (int i = 0; i < 3; ++i)
        wVecs.push_back(std::vector<double>(N_tube, 0));

    w.resize(3);

    for (int i = 0; i < w.size(); ++i)
        w[i] = &wVecs[i][0];




    // BOWING STUFF:
    fN_var = fN;
    vB_var = vB;
    //fN_var = 0.0;
    //vB_var = 0.0;
    fS_var = fN_var * mu_S;
    fC_var = fN_var * mu_C;
    z_ba_var = 0.7 * (mu_C * fN_var) / s0;
    s3_var = s3_fac * fN;

    //// What are these for ?  
    x_inp_var = bp; // bowing pos in percentage
    y_inp_var = 0.5;

    //x_out_var = 0.7;
    //y_out_var = 0.3;



    // INITIALIZE PARAMS FOR NEWTON-RAPHSON

    z_prev = 0;
    r_last = vB_var;
    r_prev = vB_var;
    a_last = r_last;
    a_prev = r_prev;

    vrel_last = -vB_var;
    z_last = 0;

    F_last = 0;


}

void MainComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    // Your audio-processing code goes here!

    // For more details, see the help for AudioProcessor::getNextAudioBlock()

    // Right now we are not producing any data, in which case we need to clear the buffer
    // (to prevent the output of random noise)
    bufferToFill.clearActiveBufferRegion();


    float* const outputL = bufferToFill.buffer->getWritePointer(0);
    float* const outputR = bufferToFill.buffer->getWritePointer(1);

    //Logger::getCurrentLogger()->outputDebugString("bufferSize: (" + juce::String(bufferToFill.numSamples) + ")");


    //AudioDeviceManager::AudioDeviceSetup currentAudioSetup;
    //deviceManager.getAudioDeviceSetup(currentAudioSetup);
    //Logger::getCurrentLogger()->outputDebugString("sample Rate in getNextAudioBlock: (" + juce::String(currentAudioSetup.sampleRate) + ")");
    //Logger::getCurrentLogger()->outputDebugString("buffer size in getNextAudioBlock: (" + juce::String(currentAudioSetup.bufferSize) + ")");



    for (int channel = 0; channel < bufferToFill.buffer->getNumChannels(); ++channel)
    {


        for (int i = 0; i < bufferToFill.numSamples; i++)
        {

            //Logger::getCurrentLogger()->outputDebugString("noSamples: (" + juce::String(bufferToFill.numSamples) + ")");
            //Logger::getCurrentLogger()->outputDebugString("i: (" + juce::String(i) + ")");

            if (channel == 0) //// this condition causes the output to explode sometimes
            {

                w_rnd = -1 + 2 * ((double)rand() / (RAND_MAX)); // random number between 0 and 1


                // Calculate bow distributions. linear interpolation in 2d.. 
                //double x_bow_ratio = x_bow_membrane / L_membrane; // can change this to x_inp_var 

                double x_bow_ratio = x_inp_var; // x_inp_var is already in percentage ! be careful between perc and absolute vals for positions..
                double y_bow_ratio = y_bow_membrane / L_membrane; 

                int l_inp = floor(x_bow_ratio * N_membrane - 1); // -1 for alignment with Matlab
                int m_inp = floor(y_bow_ratio * N_membrane - 1); // -1 for alignment with Matlab
                int l_inp_plus1 = l_inp + 1;
                int m_inp_plus1 = m_inp + 1;
                double alpha_x_inp = x_bow_ratio * N_membrane - 1 - l_inp;
                double alpha_y_inp = y_bow_ratio * N_membrane - 1 - m_inp;

                I_B[l_inp + m_inp * N_membrane] = (1 - alpha_x_inp) * (1 - alpha_y_inp); // l_inp*N+m_inp is equivalent to array[l_inp][m_inp] if it would be 2D
                I_B[l_inp_plus1 + m_inp * N_membrane] = alpha_x_inp * (1 - alpha_y_inp);
                I_B[l_inp + m_inp_plus1 * N_membrane] = (1 - alpha_x_inp) * alpha_y_inp;
                I_B[l_inp_plus1 + m_inp_plus1 * N_membrane] = alpha_x_inp * alpha_y_inp;

                for (int i = 0; i < I_B.size(); ++i)
                {
                    J_B[i] = I_B[i] * (1 / (h_membrane * h_membrane)); // speed up: keep divisions out of loop ! 
                }


                // Stuff for the update equations..
                double I_T_J_T = 0;
                double I_T_w = 0;
                double I_T_wPrev = 0;
                double I_T_dxx_w = 0;
                double I_T_dxx_wPrev = 0;
                int idx_p1;
                int idx_m1;

                for (int idx = 1; idx < N_tube - 1; ++idx)
                {
                    idx_p1 = idx + 1;
                    idx_m1 = idx - 1;

                    I_T_J_T = I_T_J_T + I_T[idx] * J_T[idx];
                    I_T_w = I_T_w + I_T[idx] * w[1][idx];
                    I_T_wPrev = I_T_wPrev + I_T[idx] * w[2][idx];
                    I_T_dxx_w = I_T_dxx_w + I_T[idx] * (w[1][idx_p1] - 2. * w[1][idx] + w[1][idx_m1]) * (1 / (h_tube * h_tube));
                    I_T_dxx_wPrev = I_T_dxx_wPrev + I_T[idx] * (w[2][idx_p1] - 2. * w[2][idx] + w[2][idx_m1]) * (1 / (h_tube * h_tube));
                }




                double I_M_J_M = 0;
                double I_B_J_B = 0;
                double I_M_J_B = 0;

                double I_M_u = 0;
                double I_M_uPrev = 0;
                double I_B_u = 0;
                double I_B_uPrev = 0;

                double I_M_dxx_u = 0;
                double I_M_dyy_u = 0;
                double I_M_delta_laplace_u = 0;
                double I_M_dxx_uPrev = 0;
                double I_M_dyy_uPrev = 0;
                double I_M_delta_laplace_uPrev = 0;

                double I_B_dxx_u = 0;
                double I_B_dyy_u = 0;
                double I_B_delta_laplace_u = 0;
                double I_B_dxx_uPrev = 0;
                double I_B_dyy_uPrev = 0;
                double I_B_delta_laplace_uPrev = 0;

                int idx_x;
                int idx_x_p1;
                int idx_x_p2;
                int idx_x_m1;
                int idx_x_m2;
                int idx_y;
                int idx_y_p1;
                int idx_y_p2;
                int idx_y_m1;
                int idx_y_m2;

                int idx_x_p1_y_p1;
                int idx_x_p1_y;
                int idx_x_p1_y_m1;
                int idx_x_y_p1;
                int idx_x_y;
                int idx_x_y_m1;
                int idx_x_m1_y_p1;
                int idx_x_m1_y;
                int idx_x_m1_y_m1;

                for (int iX = 2; iX < N_membrane - 2; ++iX)
                {
                    for (int iY = 2; iY < N_membrane - 2; ++iY)
                    {
                        idx_x = iY + (iX)*N_membrane;
                        idx_x_p1 = iY + (iX + 1) * N_membrane;
                        idx_x_p2 = iY + (iX + 2) * N_membrane;
                        idx_x_m1 = iY + (iX - 1) * N_membrane;
                        idx_x_m2 = iY + (iX - 2) * N_membrane;
                        idx_y = iY + (iX)*N_membrane;
                        idx_y_p1 = (iY + 1) + (iX)*N_membrane;
                        idx_y_p2 = (iY + 2) + (iX)*N_membrane;
                        idx_y_m1 = (iY - 1) + (iX)*N_membrane;
                        idx_y_m2 = (iY - 2) + (iX)*N_membrane;

                        idx_x_p1_y_p1 = (iY + 1) + (iX + 1) * N_membrane;
                        idx_x_p1_y = (iY)+(iX + 1) * N_membrane;
                        idx_x_p1_y_m1 = (iY - 1) + (iX + 1) * N_membrane;
                        idx_x_y_p1 = (iY + 1) + (iX)*N_membrane;
                        idx_x_y = (iY)+(iX)*N_membrane;
                        idx_x_y_m1 = (iY - 1) + (iX)*N_membrane;
                        idx_x_m1_y_p1 = (iY + 1) + (iX - 1) * N_membrane;
                        idx_x_m1_y = (iY)+(iX - 1) * N_membrane;
                        idx_x_m1_y_m1 = (iY - 1) + (iX - 1) * N_membrane;

                        I_M_J_M = I_M_J_M + I_M[idx_x] * J_M[idx_x];
                        I_B_J_B = I_B_J_B + I_B[idx_x] * J_B[idx_x];
                        I_M_J_B = I_M_J_B + I_M[idx_x] * J_B[idx_x];


                        I_M_u = I_M_u + u[1][idx_x] * I_M[idx_x];
                        I_M_uPrev = I_M_uPrev + u[2][idx_x] * I_M[idx_x];

                        I_B_u = I_B_u + u[1][idx_x] * I_B[idx_x];
                        I_B_uPrev = I_B_uPrev + u[2][idx_x] * I_B[idx_x];

                        I_M_dxx_u = I_M_dxx_u + I_M[idx_x] * (u[1][idx_x_p1] - 2. * u[1][idx_x] + u[1][idx_x_m1]) * (1 / (h_membrane * h_membrane));
                        I_M_dyy_u = I_M_dyy_u + I_M[idx_y] * (u[1][idx_y_p1] - 2. * u[1][idx_y] + u[1][idx_y_m1]) * (1 / (h_membrane * h_membrane));
                        I_M_delta_laplace_u = I_M_dxx_u + I_M_dyy_u;

                        I_M_dxx_uPrev = I_M_dxx_uPrev + I_M[idx_x] * (u[2][idx_x_p1] - 2. * u[2][idx_x] + u[2][idx_x_m1]) * (1 / (h_membrane * h_membrane));
                        I_M_dyy_uPrev = I_M_dyy_uPrev + I_M[idx_y] * (u[2][idx_y_p1] - 2. * u[2][idx_y] + u[2][idx_y_m1]) * (1 / (h_membrane * h_membrane));
                        I_M_delta_laplace_uPrev = I_M_dxx_uPrev + I_M_dyy_uPrev;

                        I_B_dxx_u = I_B_dxx_u + I_B[idx_x] * (u[1][idx_x_p1] - 2. * u[1][idx_x] + u[1][idx_x_m1]) * (1 / (h_membrane * h_membrane));
                        I_B_dyy_u = I_B_dyy_u + I_B[idx_y] * (u[1][idx_y_p1] - 2. * u[1][idx_y] + u[1][idx_y_m1]) * (1 / (h_membrane * h_membrane));
                        I_B_delta_laplace_u = I_B_dxx_u + I_B_dyy_u;

                        I_B_dxx_uPrev = I_B_dxx_uPrev + I_B[idx_x] * (u[2][idx_x_p1] - 2. * u[2][idx_x] + u[2][idx_x_m1]) * (1 / (h_membrane * h_membrane));
                        I_B_dyy_uPrev = I_B_dyy_uPrev + I_B[idx_y] * (u[2][idx_y_p1] - 2. * u[2][idx_y] + u[2][idx_y_m1]) * (1 / (h_membrane * h_membrane));
                        I_B_delta_laplace_uPrev = I_B_dxx_uPrev + I_B_dyy_uPrev;

                        //I_P_dxxxx_w = I_P_dxxxx_w + I_P[idx_x] * (w[1][idx_x_p2] - 4. * w[1][idx_x_p1] + 6. * w[1][idx_x] - 4. * w[1][idx_x_m1] + w[1][idx_x_m2]) * (1 / (h_plate * h_plate * h_plate * h_plate));
                        //I_P_dyyyy_w = I_P_dyyyy_w + I_P[idx_y] * (w[1][idx_y_p2] - 4. * w[1][idx_y_p1] + 6. * w[1][idx_y] - 4. * w[1][idx_y_m1] + w[1][idx_y_m2]) * (1 / (h_plate * h_plate * h_plate * h_plate));
                        //x2_I_P_dxxyy_w = x2_I_P_dxxyy_w + I_P[idx_x]
                        //    * (w[1][idx_x_p1_y_p1] - 2. * w[1][idx_x_p1_y] + w[1][idx_x_p1_y_m1] - 2. * w[1][idx_x_y_p1] + 4. * w[1][idx_x_y]
                        //        - 2. * w[1][idx_x_y_m1] + w[1][idx_x_m1_y_p1] - 2. * w[1][idx_x_m1_y] + w[1][idx_x_m1_y_m1]) * (2 / (h_plate * h_plate * h_plate * h_plate));

                        //I_P_delta_laplace_x2_w = I_P_dxxxx_w + I_P_dyyyy_w + x2_I_P_dxxyy_w;
                    }
                }

                // Can re - write q so it is clearer..

                // Terms for the Newton-Raphson:
                double q = -(1 / (fN_var * I_B_J_B)) * ((-2 / k) * (1 / k) * (I_B_u - I_B_uPrev)
                        + (2 / k) * vB
                        + 2 * sig0 * vB
                        - (c_membrane * c_membrane) * I_B_delta_laplace_u
                        + (2 * sig1) * (1 / k) * I_B_delta_laplace_u
                        + (2 * sig1) * (1 / k) * I_B_delta_laplace_uPrev);

                double b = (k * k / (1 + sig0 * k)) * (c_membrane * c_membrane * I_M_delta_laplace_u + (2 * sig1 / k) * (I_M_delta_laplace_u - I_M_delta_laplace_uPrev)) 
                    + (k * k / (1 + sig0 * k)) * (2 / (k * k)) * I_M_u - ((k * k) / (1 + sig0 * k)) * ((1 - sig0 * k) / (k * k)) * I_M_uPrev
                    - (k * k / (1 + sig0_tube * k)) * (c_tube * c_tube * I_T_dxx_w + (2 * sig1_tube / k) * (I_T_dxx_w - I_T_dxx_wPrev))
                    - (k * k / (1 + sig0_tube * k)) * (2 / (k * k)) * I_T_w + ((k * k) / (1 + sig0_tube * k)) * ((1 - sig0_tube * k) / (k * k)) * I_T_wPrev;



                s3 = s3_fac * fN_var; // fN_var is FB in Matlab code..
                q = q / -(1 / (fN_var * I_B_J_B)); // hope this is good ? !where the hell did this come from ? : /

                // Newton-Raphson Solver for bowing model. Just do it on 2 params first and include the connection force later ! 
                double eps = 1;
                int iter_check = 0;
                vrel_last = -vB;
                z_last = 0;
                // Newton-Raphson iterative scheme
                while ((eps > tol) && (fC_var > 0))
                {
                    ++iter_check;

                    if (vrel_last == 0)
                    {
                        z_ss = fS_var / s0;
                    }
                    else
                    {
                        z_ss = sign(vrel_last) / s0 * (fC_var + (fS_var - fC_var) * exp(-pow((vrel_last / vS), 2)));
                    }

                    if (sign(vrel_last) == sign(z_last))
                    {
                        if (abs(z_last) <= z_ba_var)
                        {
                            alpha_fct = 0;
                        }
                        else if ((z_ba_var < abs(z_last)) && (abs(z_last) < abs(z_ss)))
                        {
                            alpha_fct = 0.5 * (1 + sign(z_last) * sin((double_Pi * (z_last - sign(z_last) * 0.5 * (abs(z_ss) + z_ba_var)) / (abs(z_ss) - z_ba_var))));
                        }
                        else if (abs(z_last) >= abs(z_ss))
                        {
                            alpha_fct = 1;
                        }
                    }
                    else
                    {
                        alpha_fct = 0;
                    }


                    r_last = vrel_last * (1 - alpha_fct * z_last / z_ss);
                    a_last = ((2 / k) * (z_last - z_prev) - a_prev);

                    g1 = I_B_J_B * ((s0 * z_last + s1 * r_last + s2 * vrel_last + s3_var * w_rnd) / (rho_membrane * H_membrane)) + (2 / k + 2 * sig0) * vrel_last + q;
                    g2 = r_last - ((2 / k) * (z_last - z_prev) - a_prev);

                    if (sign(vrel_last) >= 0)
                    {
                        z_ss_d_vrel = -2 * vrel_last * (-fC_var + fS_var) * exp(-(vrel_last * vrel_last) / (vS * vS)) / (s0 * (vS * vS));
                    }
                    else
                    {
                        z_ss_d_vrel = 2 * vrel_last * (-fC_var + fS_var) * exp(-(vrel_last * vrel_last) / (vS * vS)) / (s0 * (vS * vS));
                    }

                    if ((z_ba_var < abs(z_last)) && (abs(z_last) < abs(z_ss)))
                    {
                        if (sign(z_last) >= 0)
                        {
                            alpha_fct_d_vrel = 0.5 * (-0.5 * double_Pi * (z_ss * z_ss_d_vrel) * sign(z_ss) / ((-z_ba_var + abs(z_ss)) * z_ss) - double_Pi * (z_ss * z_ss_d_vrel) * (-0.5 * z_ba_var + z_last - 0.5 * abs(z_ss)) * sign(z_ss) / (pow((-z_ba_var + abs(z_ss)), 2) * z_ss)) * cos(double_Pi * (-0.5 * z_ba_var + z_last - 0.5 * abs(z_ss)) / (-z_ba_var + abs(z_ss)));
                            alpha_fct_d_z = 0.5 * double_Pi * cos(double_Pi * (-0.5 * z_ba_var + z_last - 0.5 * abs(z_ss)) / (-z_ba_var + abs(z_ss))) / (-z_ba_var + abs(z_ss));
                        }
                        else
                        {
                            alpha_fct_d_vrel = -0.5 * (0.5 * double_Pi * (z_ss * z_ss_d_vrel) * sign(z_ss) / ((-z_ba_var + abs(z_ss)) * z_ss) - double_Pi * (z_ss * z_ss_d_vrel) * (0.5 * z_ba_var + z_last + 0.5 * abs(z_ss)) * sign(z_ss) / (pow((-z_ba_var + abs(z_ss)), 2) * z_ss)) * cos(double_Pi * (0.5 * z_ba_var + z_last + 0.5 * abs(z_ss)) / (-z_ba_var + abs(z_ss)));
                            alpha_fct_d_z = -0.5 * double_Pi * cos(double_Pi * (0.5 * z_ba_var + z_last + 0.5 * abs(z_ss)) / (-z_ba_var + abs(z_ss))) / (-z_ba_var + abs(z_ss));
                        }
                    }
                    else
                    {
                        alpha_fct_d_vrel = 0;
                        alpha_fct_d_z = 0;
                    }

                    r_last_d_vrel = vrel_last * (z_last * alpha_fct * z_ss_d_vrel / (z_ss * z_ss) - z_last * alpha_fct_d_vrel / z_ss) - z_last * alpha_fct / z_ss + 1;
                    r_last_d_z = vrel_last * (-z_last * alpha_fct_d_z / z_ss - alpha_fct / z_ss);

                    g1_d_vrel = 2 * sig0 + 2 / k + I_B_J_B * (s1 * r_last_d_vrel + s2) / (rho_membrane * H_membrane);
                    g1_d_z = I_B_J_B * (s0 + s1 * r_last_d_z) / (rho_membrane * H_membrane);

                    g2_d_vrel = r_last_d_vrel;
                    g2_d_z = r_last_d_z - 2 / k;

                    fac_1_over_det = 1 / (g1_d_vrel * g2_d_z - g1_d_z * g2_d_vrel);

                    A_mat_new = fac_1_over_det * g2_d_z;
                    B_mat_new = fac_1_over_det * (-g1_d_z);
                    C_mat_new = fac_1_over_det * (-g2_d_vrel);
                    D_mat_new = fac_1_over_det * (g1_d_vrel);

                    vrel = vrel_last - (A_mat_new * g1 + B_mat_new * g2);
                    z = z_last - (C_mat_new * g1 + D_mat_new * g2);

                    eps = sqrt((z - z_last) * (z - z_last) + (vrel - vrel_last) * (vrel - vrel_last)); // same as norm(theta - theta_last);

                    z_last = z;
                    vrel_last = vrel;

                    if (iter_check == 99)
                        break;
                }

                if (fC_var > 0)
                {
                    if (vrel == 0)
                    {
                        z_ss = fS_var / s0;
                    }
                    else
                    {
                        //z_ss = sign(vrel) / s0 * (fC_var + (fS_var - fC_var) * exp(-pow((vrel / vS), 2)));
                        z_ss = sign(vrel) / s0 * (fC_var + (fS_var - fC_var) * exp(-(vrel / vS) * (vrel / vS)));
                    }

                    double ana = sign(vrel);

                    if (sign(vrel) == sign(z))
                    {
                        if (abs(z) <= z_ba_var)
                        {
                            alpha_fct = 0;
                        }
                        else if ((z_ba_var < abs(z)) && (abs(z) < abs(z_ss)))
                        {
                            alpha_fct = 0.5 * (1 + sign(z) * sin((double_Pi * (z - sign(z) * 0.5 * (abs(z_ss) + z_ba_var)) / (abs(z_ss) - z_ba_var))));
                        }
                        else if (abs(z) >= abs(z_ss))
                        {
                            alpha_fct = 1;
                        }
                    }
                    else
                    {
                        alpha_fct = 0;
                    }

                    r = vrel * (1 - alpha_fct * z / z_ss);
                    a = (2 / k) * (z - z_prev) - a_prev;


                    Fbow = (s0 * z + s1 * r + s2 * vrel + s3_var * w_rnd);


                    z_prev = z;
                    a_prev = a;
                }
                else
                {
                    Fbow = 0;

                    // reset params for NR
                    z_prev = 0;
                    r_last = vB_var;
                    r_prev = vB_var;
                    a_last = r_last;
                    a_prev = r_prev;

                    vrel_last = -vB_var;
                    z_last = 0;

                    F_last = 0;
                }

                // I'm not sure this is correct..
                F = (-b + (k * k / (1 + sig0 * k)) * I_M_J_B * ((s0 * z_last + s1 * r_last + s2 * vrel_last + s3 * w_rnd) / (rho_membrane * H_membrane))) / (k * k * I_T_J_T / ((1 + sig0 * k) * rho_tube * A_tube) + k * k * I_M_J_M / ((1 + sig0_tube * k) * rho_tube * A_tube));


                // Update equations:

                for (int iLoc = 0; iLoc < locsDo.size(); ++iLoc) 
                {
                    int iX = locsDo[iLoc][0];
                    int iY = locsDo[iLoc][1];

                    idx_x = iY + (iX)*N_membrane;
                    idx_x_p1 = iY + (iX + 1) * N_membrane;
                    idx_x_p2 = iY + (iX + 2) * N_membrane;
                    idx_x_m1 = iY + (iX - 1) * N_membrane;
                    idx_x_m2 = iY + (iX - 2) * N_membrane;
                    idx_y = iY + (iX)*N_membrane;
                    idx_y_p1 = (iY + 1) + (iX)*N_membrane;
                    idx_y_p2 = (iY + 2) + (iX)*N_membrane;
                    idx_y_m1 = (iY - 1) + (iX)*N_membrane;
                    idx_y_m2 = (iY - 2) + (iX)*N_membrane;

                    idx_x_p1_y_p1 = (iY + 1) + (iX + 1) * N_membrane;
                    idx_x_p1_y = (iY)+(iX + 1) * N_membrane;
                    idx_x_p1_y_m1 = (iY - 1) + (iX + 1) * N_membrane;
                    idx_x_y_p1 = (iY + 1) + (iX)*N_membrane;
                    idx_x_y = (iY)+(iX)*N_membrane;
                    idx_x_y_m1 = (iY - 1) + (iX)*N_membrane;
                    idx_x_m1_y_p1 = (iY + 1) + (iX - 1) * N_membrane;
                    idx_x_m1_y = (iY)+(iX - 1) * N_membrane;
                    idx_x_m1_y_m1 = (iY - 1) + (iX - 1) * N_membrane;


                    u[0][idx_x_y] = (k * k / (1 + sig0 * k)) * ((c_membrane * c_membrane / (h_membrane * h_membrane)) * (u[1][idx_x_p1] + u[1][idx_x_m1] + u[1][idx_y_p1] + u[1][idx_y_m1] - 4 * u[1][idx_x_y])
                        + (2 * sig1 / k) * (1 / (h_membrane * h_membrane)) * (u[1][idx_x_p1] + u[1][idx_x_m1] + u[1][idx_y_p1] + u[1][idx_y_m1] - 4 * u[1][idx_x_y] - (u[2][idx_x_p1] + u[2][idx_x_m1] + u[2][idx_y_p1] + u[2][idx_y_m1] - 4 * u[2][idx_x_y]))
                        - J_B[idx_x_y] * Fbow / (rho_membrane * H_membrane) + J_M[idx_x_y] * F / (rho_membrane * H_membrane)
                        + (2 / (k * k)) * u[1][idx_x_y] - (1 - sig0 * k) * u[2][idx_x_y] / (k * k));

                }


                //// Calculate on the whole grid! 
                //for (int iX = 1; iX < N_membrane - 1; ++iX)
                //{
                //    for (int iY = 1; iY < N_membrane - 1; ++iY)
                //    {
                //        idx_x = iY + (iX)*N_membrane;
                //        idx_x_p1 = iY + (iX + 1) * N_membrane;
                //        idx_x_p2 = iY + (iX + 2) * N_membrane;
                //        idx_x_m1 = iY + (iX - 1) * N_membrane;
                //        idx_x_m2 = iY + (iX - 2) * N_membrane;
                //        idx_y = iY + (iX)*N_membrane;
                //        idx_y_p1 = (iY + 1) + (iX)*N_membrane;
                //        idx_y_p2 = (iY + 2) + (iX)*N_membrane;
                //        idx_y_m1 = (iY - 1) + (iX)*N_membrane;
                //        idx_y_m2 = (iY - 2) + (iX)*N_membrane;

                //        idx_x_p1_y_p1 = (iY + 1) + (iX + 1) * N_membrane;
                //        idx_x_p1_y = (iY)+(iX + 1) * N_membrane;
                //        idx_x_p1_y_m1 = (iY - 1) + (iX + 1) * N_membrane;
                //        idx_x_y_p1 = (iY + 1) + (iX)*N_membrane;
                //        idx_x_y = (iY)+(iX)*N_membrane;
                //        idx_x_y_m1 = (iY - 1) + (iX)*N_membrane;
                //        idx_x_m1_y_p1 = (iY + 1) + (iX - 1) * N_membrane;
                //        idx_x_m1_y = (iY)+(iX - 1) * N_membrane;
                //        idx_x_m1_y_m1 = (iY - 1) + (iX - 1) * N_membrane;


                //        u[0][idx_x_y] = (k * k / (1 + sig0 * k)) * ((c_membrane * c_membrane / (h_membrane * h_membrane)) * (u[1][idx_x_p1] + u[1][idx_x_m1] + u[1][idx_y_p1] + u[1][idx_y_m1] - 4 * u[1][idx_x_y])
                //                      + (2 * sig1 / k) * (1 / (h_membrane * h_membrane)) * (u[1][idx_x_p1] + u[1][idx_x_m1] + u[1][idx_y_p1] + u[1][idx_y_m1] - 4 * u[1][idx_x_y] - (u[2][idx_x_p1] + u[2][idx_x_m1] + u[2][idx_y_p1] + u[2][idx_y_m1] - 4 * u[2][idx_x_y])) 
                //                      - J_B[idx_x_y] * Fbow / (rho_membrane * H_membrane) + J_M[idx_x_y] * F / (rho_membrane * H_membrane)
                //                      + (2 / (k * k)) * u[1][idx_x_y] - (1 - sig0 * k) * u[2][idx_x_y] / (k * k));


                //    }
                //}


                int idx_p2;
                int idx_m2;
                for (int idx = 1; idx < N_tube - 1; ++idx)
                {
                    idx_p1 = idx + 1;
                    idx_p2 = idx + 2;
                    idx_m1 = idx - 1;
                    idx_m2 = idx - 2;

                    // maybe you can reused the definitions from above I_S_etc

                    w[0][idx]  = (k * k) * ( (c_tube * c_tube / (h_tube * h_tube)) * (w[1][idx_p1] - 2 * w[1][idx] + w[1][idx_m1])
                        - J_T[idx] * F / (rho_tube * A_tube) 
                        + (2 / (k * k)) * w[1][idx] - w[2][idx] / (k * k) );
                }

                w[0][0] = (k * k) * ( (c_tube * c_tube / (h_tube * h_tube)) * (2 * w[1][1] - 2 * w[1][0])
                    - J_T[0] * F / (rho_tube * A_tube)
                    + (2 / (k * k)) * w[1][0] - w[2][0] / (k * k));

                w[0][(int)N_tube - 1] = (w[1][(int)N_tube - 1] * F2
                    + w[2][(int)N_tube - 1] * (F1 * (alpha1 * h_tube / k - alpha2 * h_tube) - 1)
                    + w[1][(int)N_tube - 2] * (F1 + F3)) / (F1 * (alpha1 * h_tube / k
                    + alpha2 * h_tube) + 1) 
                    - (k * k) * J_T[(int)N_tube - 1] *F / (rho_tube * A_tube);



                preOutput = w[0][(int)N_tube - 1]; // take the output at the end of the tube..
                preOutput = preOutput * gain;


                uTmp = u[2];
                u[2] = u[1];
                u[1] = u[0];
                u[0] = uTmp;

                wTmp = w[2];
                w[2] = w[1];
                w[1] = w[0];
                w[0] = wTmp;


                //uVecs

                // Pointers pointers pointers !!! 
                // Give you more control on what your data does and how to point to stuff instead of copying stuff.

                //if (i == 0)
                //{
                //    Logger::getCurrentLogger()->outputDebugString("preOutput: (" + String(preOutput) + ")");
                //}

                //preOutput = gain * sineWave;
                if (abs(preOutput) > 1)
                {
                    /*std::cout << "Output is too loud!" << std::endl;*/
                    //Logger::getCurrentLogger()->outputDebugString("Output is too loud!");
                }
                else {
                    //output_interim[i] = preOutput;
                    outputL[i] = preOutput;
                    ////outputR[i] = gain * uNext[floor(6 * N / 8)];
                    outputR[i] = outputL[i];
                }

                //uVecs

                //Logger::getCurrentLogger()->outputDebugString(String(t));
                ++t;
                int anabanana = 5;


            }
        }
    }


}

void MainComponent::releaseResources()
{
    // This will be called when the audio device stops, or when it is being
    // restarted due to a setting change.

    // For more details, see the help for AudioProcessor::releaseResources()
}

//==============================================================================
void MainComponent::paint (juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (juce::ResizableWindow::backgroundColourId));

    // You can add your drawing code here!
}

void MainComponent::resized()
{
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.


    auto sliderLeft = 120;
    frequencySlider.setBounds(sliderLeft, 20, getWidth() - sliderLeft - 10, 20);

    sig_0_Slider.setBounds(sliderLeft, 50, getWidth() - sliderLeft - 10, 20);

    sig_1_Slider.setBounds(sliderLeft, 80, getWidth() - sliderLeft - 10, 20);

}
