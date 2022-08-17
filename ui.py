import numpy as np
import streamlit as st
from functions import run_computation
import matplotlib.pyplot as plt
from functions import run_computation
def run_(r_outer, r_inner, r, length, a_min, a_max, d_max, d_min, freq, phi, C):
    # clear st.session_state.figure
    initialize()
    I = np.eye(3,3)
    phi = np.radians(phi)
    # TODO run computations
    r = r_outer if r == "Outer Radius" else r_inner
    Rez, N0_crit, CRIT, sig_maxtt  = run_computation(r_outer, r_inner, r, length, a_min, a_max, d_max, d_min, freq, phi, c_1=C)
    values = []
    values.append(N0_crit)
    values.append(CRIT)
    values.append(sig_maxtt)
    fig = []
    title = []
    plt.rcParams['font.size'] = '7'
    fig1 = plt.figure()
    fig2 = plt.figure()
    fig3 = plt.figure()
    fig4 = plt.figure()
    fig5 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(Rez[:,0], Rez[:,6], color='r', label="sigma_zz")
    ax1.plot(Rez[:,0], Rez[:,7], color='b', label="sigma_thetaz")
    title.append("Evolution of axial and shear Cauchy stresses")
    ax1.legend()  
    fig.append(fig1)  
    ax2 = fig2.add_subplot(111)
    ax2.plot(Rez[:,0], Rez[:,10], color="r", label="sigma_n")
    ax2.plot(Rez[:,0], Rez[:,11], color="b", label="tau_s")
    title.append("Evolution of normal and tangential components of traction vector")
    ax2.legend()
    fig.append(fig2)
    ax3 = fig3.add_subplot(111)
    ax3.plot(Rez[:,0], Rez[:,13], color="r", label="sigma_1")
    ax3.plot(Rez[:,0], Rez[:,14], color="b", label="sigma_2")
    ax3.plot(Rez[:,0], Rez[:,15], color="g", label="sigma_3")
    title.append("Evolution of the principal Cauchy stresses")
    ax3.legend()
    fig.append(fig3)
    ax4 = fig4.add_subplot(111)
    ax4.plot(Rez[:,0], Rez[:,16], color="r", label="Sigma_1")
    ax4.plot(Rez[:,0], Rez[:,17], color="b", label="Sigma_2")
    ax4.plot(Rez[:,0], Rez[:,18], color="g", label="Sigma_3")
    title.append("Evolution of the principal configurational stresses")
    ax4.legend()
    fig.append(fig4)
    ax5 = fig5.add_subplot(111)
    ax5.plot(Rez[:,0], Rez[:,19], color="r", label="sigma_ZZ")
    ax5.plot(Rez[:,0], Rez[:,20], color="b", label="sigma_thetaZ")
    ax5.plot(Rez[:,0], Rez[:,21], color="g", label="sigma_RR")
    title.append("Evolution of axial, shear and radial configurational stresses")
    ax5.legend()
    fig.append(fig5)
    plt.subplots_adjust(hspace=0.7, wspace=1.8)
    st.session_state.figure = fig
    st.session_state.title = title
    st.session_state.v = values

def generate_sidebar(sidebar:st.sidebar):
    sidebar.title("Inputs")
    sidebar.markdown("#")
    sidebar.markdown("#")
    geometrics = sidebar.expander("Cylinder Geometrics")
    r_outer = geometrics.number_input("Outer Radius - millimeters", value=43.18)
    r_inner = geometrics.number_input("Inner Radius - millimeters", value=38.1)
    r = geometrics.selectbox("R", options=["Outer Radius", "Inner Radius"])
    length = geometrics.number_input("Cylinder Length - millimeters", value=6.35)
    sidebar.markdown("#")
    comp_conf = sidebar.expander("Computation Configuration")
    a_min = comp_conf.number_input("Minimum Twist Angle - Degrees", value=0)
    a_max = comp_conf.number_input("Maximum Twist Angle - Degrees", value=6)
    d_min = comp_conf.number_input("Minimum Displacement - millimeters", value=0)
    d_max = comp_conf.number_input("Maximum Displacement - millimeters", value=0.5*length)
    freq = comp_conf.number_input("Frequency - Hz", value=10)
    phi = comp_conf.number_input("Phase Angle - Degrees", value=0)
    sidebar.markdown("#")
    materials = sidebar.expander("Cylinder Material")
    C = materials.number_input("Material - MPa", value=1.5)
    sidebar.markdown("#")
    sidebar.markdown("#")
    run = sidebar.button(
        label="Run Computation", 
        key="run",
        on_click=run_, 
        args=(r_outer, r_inner, r, length, a_min, a_max, d_max, d_min, freq, phi, C))
    sidebar.write(" ")

def initialize():
    st.session_state.v = []
    st.session_state.figure = []
    st.session_state.title = []

def main():
    sidebar = st.sidebar
    generate_sidebar(sidebar)
    default = "None"
    st.title("Traction Tortion Application")
    st.header("Results")
    st.markdown("#")
    st.markdown("#")
    if "figure" in st.session_state:
        if st.session_state.figure:
            title = st.session_state.title
            for i, fig in enumerate(st.session_state.figure):
                st.markdown(f"<h4 style = 'text-align: center'>{title[i]}<h2>", unsafe_allow_html=True)
                st.pyplot(fig)
        values = st.session_state.v
        col_1, col_2, col_3 = st.columns(3)
        col_1.write("N0 crit: "+ str(values[0] if values else default))
        col_2.write("CRIT: "+ str(values[1] if values else default))
        col_3.write("sig_maxtt: "+ str(values[2] if values else default))
    else:
        st.write("Run a computation for the graph to be displayed.")


main()