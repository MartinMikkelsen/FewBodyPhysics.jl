module Constants

using Unitful
export m_p, m_n, m_pi0, m_pi, ħc, m_bare, μ

const uMeVc2 = u"MeV*1/c^2"

const m_p   = 938.27uMeVc2       # Proton mass
const m_n   = 939.57uMeVc2       # Neutron mass
const m_pi0 = 134.98uMeVc2       # Neutral pion
const m_pi  = 139.57uMeVc2       # Charged pion

const ħc    = 197.3269804u"MeV*fm"  # ħ·c in MeV·fm

# Derived masses
const m_bare = (m_p + m_n) / 2        # Average nucleon mass
const μ = (m_bare * m_pi0) / (m_bare + m_pi0)  # Reduced mass of pion-nucleon system

end
