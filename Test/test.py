import torch
import pyro
import pyro.distributions as dist
from pyro.infer import SVI, Trace_ELBO
from pyro.optim import Adam
from pyro.contrib.hmm import DiscreteHMM

# Parameters
K = 2  # Number of states
T = 52.18  # Period
M = 2  # Number of mixture components

def compute_Pij(t, beta_ij):
    # Convert t to tensor if it's not already
    t_tensor = torch.tensor(float(t), dtype=torch.float32)
    
    # Only use the first 2 harmonics (indices 0, 1, 2, 3, 4)
    # beta_ij[0] is the constant term
    # beta_ij[1] is cos(2πt/T), beta_ij[2] is sin(2πt/T)
    # beta_ij[3] is cos(4πt/T), beta_ij[4] is sin(4πt/T)
    return (beta_ij[0] + 
            beta_ij[1] * torch.cos(2 * torch.pi * t_tensor / T) +
            beta_ij[2] * torch.sin(2 * torch.pi * t_tensor / T) +
            beta_ij[3] * torch.cos(4 * torch.pi * t_tensor / T) +
            beta_ij[4] * torch.sin(4 * torch.pi * t_tensor / T))

def transition_matrix(t, beta):
    # Use softmax for proper probability distribution
    P = torch.zeros((K, K))
    for i in range(K):
        logits = torch.zeros(K)
        for j in range(K-1):
            logits[j+1] = compute_Pij(t, beta[i, j])
        P[i] = torch.softmax(logits, dim=0)
    return P

def model(data_mortality, data_temperature, exposures):
    # Sample global parameters with priors
    beta_prior = torch.randn(K, K-1, 5) * 0.1
    beta = pyro.sample('beta', dist.Normal(beta_prior, 0.1).to_event(3))
    
    lambda_km = pyro.sample('lambda_km', dist.Normal(torch.zeros(K, M), 1.0).to_event(2))
    alpha_k = pyro.sample('alpha_k', dist.Normal(torch.zeros(K), 0.1).to_event(1))
    
    # gamma_k should have 4 components for the 2 harmonics
    # [cos(2πt/T), sin(2πt/T), cos(4πt/T), sin(4πt/T)]
    gamma_k = pyro.sample('gamma_k', dist.Normal(torch.zeros(K, 4), 0.1).to_event(2))
    
    # Mixture weights (using softmax for proper probabilities)
    p_km_logits = pyro.sample('p_km_logits', dist.Normal(torch.zeros(K, M), 1.0).to_event(2))
    p_km = torch.softmax(p_km_logits, dim=1)
    
    # Parameters for Normal component distributions
    m_km = pyro.sample('m_km', dist.Normal(torch.zeros(K, M), 5.0).to_event(2))
    sigma_km = pyro.sample('sigma_km', dist.LogNormal(torch.zeros(K, M), 0.5).to_event(2))

    # Initial state distribution
    state_probs = torch.ones(K) / K
    
    # Transition matrix for all time steps
    transition_matrices = torch.stack([transition_matrix(t, beta) for t in range(len(data_mortality))])
    
    # Emission distributions for mortality and temperature
    mortality_dists = []
    temperature_dists = []
    for t in range(len(data_mortality)):
        t_tensor = torch.tensor(float(t), dtype=torch.float32)
        
        # Calculate mortality rate components
        T1_kt = alpha_k  # trend component
        S1_kt = (gamma_k[:, 0] * torch.cos(2 * torch.pi * t_tensor / T) +
                 gamma_k[:, 1] * torch.sin(2 * torch.pi * t_tensor / T) +
                 gamma_k[:, 2] * torch.cos(4 * torch.pi * t_tensor / T) +
                 gamma_k[:, 3] * torch.sin(4 * torch.pi * t_tensor / T))
        
        # Mortality rates for each mixture component
        log_rates = lambda_km + T1_kt.unsqueeze(1) + S1_kt.unsqueeze(1)
        rates = torch.exp(log_rates) * exposures[t]
        
        # Mixture of Poissons for mortality
        mortality_dist = dist.MixtureSameFamily(
            dist.Categorical(p_km),
            dist.Poisson(rates)
        )
        mortality_dists.append(mortality_dist)
        
        # Temperature model - mixture of normals
        T2_kt = alpha_k  # reusing alpha for simplicity
        S2_kt = (gamma_k[:, 0] * torch.cos(2 * torch.pi * t_tensor / T) +
                 gamma_k[:, 1] * torch.sin(2 * torch.pi * t_tensor / T) +
                 gamma_k[:, 2] * torch.cos(4 * torch.pi * t_tensor / T) +
                 gamma_k[:, 3] * torch.sin(4 * torch.pi * t_tensor / T))
        
        means = m_km + T2_kt.unsqueeze(1) + S2_kt.unsqueeze(1)
        temperature_dist = dist.MixtureSameFamily(
            dist.Categorical(p_km),
            dist.Normal(means, sigma_km)
        )
        temperature_dists.append(temperature_dist)
    
    # Combine mortality and temperature distributions
    emission_dists = [dist.Independent(dist.JointDistribution([mortality_dists[t], temperature_dists[t]]), 1) for t in range(len(data_mortality))]
    
    # Define the HMM
    hmm = DiscreteHMM(state_probs, transition_matrices, emission_dists)
    
    # Observe data
    pyro.sample("obs", hmm, obs=torch.stack([data_mortality, data_temperature], dim=-1))

# Define guide (AutoDiagonalNormal would be better in practice)
def guide(data_mortality, data_temperature, exposures):
    # Variational parameters for global parameters
    beta_loc = pyro.param('beta_loc', torch.zeros(K, K-1, 5))
    beta_scale = pyro.param('beta_scale', torch.ones(K, K-1, 5) * 0.1, 
                           constraint=dist.constraints.positive)
    beta = pyro.sample('beta', dist.Normal(beta_loc, beta_scale).to_event(3))
    
    lambda_loc = pyro.param('lambda_loc', torch.zeros(K, M))
    lambda_scale = pyro.param('lambda_scale', torch.ones(K, M) * 0.1,
                             constraint=dist.constraints.positive)
    lambda_km = pyro.sample('lambda_km', dist.Normal(lambda_loc, lambda_scale).to_event(2))
    
    alpha_loc = pyro.param('alpha_loc', torch.zeros(K))
    alpha_scale = pyro.param('alpha_scale', torch.ones(K) * 0.1,
                            constraint=dist.constraints.positive)
    alpha_k = pyro.sample('alpha_k', dist.Normal(alpha_loc, alpha_scale).to_event(1))
    
    gamma_loc = pyro.param('gamma_loc', torch.zeros(K, 4))
    gamma_scale = pyro.param('gamma_scale', torch.ones(K, 4) * 0.1,
                            constraint=dist.constraints.positive)
    gamma_k = pyro.sample('gamma_k', dist.Normal(gamma_loc, gamma_scale).to_event(2))
    
    p_logits_loc = pyro.param('p_logits_loc', torch.zeros(K, M))
    p_logits_scale = pyro.param('p_logits_scale', torch.ones(K, M) * 0.1,
                               constraint=dist.constraints.positive)
    p_km_logits = pyro.sample('p_km_logits', dist.Normal(p_logits_loc, p_logits_scale).to_event(2))
    
    m_loc = pyro.param('m_loc', torch.zeros(K, M))
    m_scale = pyro.param('m_scale', torch.ones(K, M),
                        constraint=dist.constraints.positive)
    m_km = pyro.sample('m_km', dist.Normal(m_loc, m_scale).to_event(2))
    
    sigma_loc = pyro.param('sigma_loc', torch.zeros(K, M))
    sigma_scale = pyro.param('sigma_scale', torch.ones(K, M) * 0.1,
                            constraint=dist.constraints.positive)
    sigma_km = pyro.sample('sigma_km', dist.LogNormal(sigma_loc, sigma_scale).to_event(2))

# Example data
data_mortality = torch.tensor([5., 10., 7., 3., 12.])
data_temperature = torch.tensor([15.3, 16.1, 14.8, 13.5, 17.2])
exposures = torch.tensor([100., 120., 110., 90., 130.])

# Initialize Pyro
pyro.clear_param_store()

# Setup optimizer and SVI
optimizer = Adam({"lr": 0.01})
svi = SVI(model, guide, optimizer, loss=Trace_ELBO())

# Training loop
n_steps = 500
for step in range(n_steps):
    loss = svi.step(data_mortality, data_temperature, exposures)
    if step % 50 == 0:
        print(f'Step {step} - Loss: {loss}')

# Get estimated parameters
for name, value in pyro.get_param_store().items():
    print(f"{name}: {value}")