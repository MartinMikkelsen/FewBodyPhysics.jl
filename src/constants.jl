export m_pi, m_pi0, m_p, m_n, μ, ħc, m_bare

const m_p   = 938.27       # Proton mass in MeV/c^2
const m_n   = 939.57       # Neutron mass in MeV/c^2
const m_pi0 = 134.98       # Neutral pion mass (π⁰) in MeV/c^2
const m_pi  = 139.57       # Charged pion mass (π⁺/π⁻) in MeV/c^2
const ħc    = 197.3269804  # Planck's constant times speed of light in MeV·fm

const m_bare = (m_p + m_n) / 2    # Average nucleon mass in MeV/c^2

"""
    μ = (m_bare * m_pi0) / (m_bare + m_pi0)

Reduced mass of the nucleon-pion system in MeV/c^2.
"""
const μ = (m_bare * m_pi0) / (m_bare + m_pi0)
