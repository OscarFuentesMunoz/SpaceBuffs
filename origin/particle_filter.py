import numpy as np
from numpy.random import default_rng
import pykep as pk
import time
import multiprocessing as mp
from scipy.stats import multivariate_normal
from origin.propagation import Propagator
from origin.const import GMe
from origin.util import normalize_diff



# Implement Bootstrap PF where the importance pdf is simply chosen to be the prior p(x_{k} | x_{k-1}).
# Arulampalam, M. Sanjeev, et al. "A tutorial on particle filters for online nonlinear/non-Gaussian Bayesian tracking." IEEE Transactions on signal processing 50.2 (2002): 174-188.
# The algorithm resamples the particles based on "resample-move" algorithm using random walk Metropolis-Hastings sampling.
class PF_Bootstrap_OE():
    """
    PF_Bootstrap_OE assumes the state is [a,e,i,W,w,E,cram] and measurement is [a,e,i,W,w,E].
    """
    def __init__(self, Q=np.zeros((7,7)), R0=np.ones((6,6)), R1=np.ones((6,6)), N=50, switch_R_after=2) -> None:
        self.Q = Q
        self.R0 = R0
        self.R1 = R1
        self.R = R0
        self.N = N
        self.switch_R_after = switch_R_after
        self.rng = default_rng(1)
        self.log_cram_range = np.array((-0.5, 1.8))
    
    def __call__(self, time_obs, obs_his, oe0, cov_oe0, is_parallel=False):
        # Initialize 
        prop = Propagator()
        weight = np.ones(self.N) / self.N
        
        # Initialize particles
        samples = np.zeros((len(time_obs), self.N, 7))
        # Add Gaussian noise to the initial state
        samples[0, :, 0:6] = self.rng.multivariate_normal(oe0, cov_oe0, (self.N, ))
        # As for Cr(A/m), sample log(cram) from a uniform distribution
        log_cram_samples = (self.log_cram_range[1] - self.log_cram_range[0]) * self.rng.random(size=self.N) + self.log_cram_range[0]
        samples[0, :, 6] = 10**log_cram_samples

        # Output array
        mean_out = np.zeros((len(time_obs), 7))
        cov_diag_out = np.zeros((len(time_obs), 7))
        mean, cov_diag = self._calc_mean_cov(samples[0, :, :])
        mean_out[0, :] = mean
        cov_diag_out[0, :] = cov_diag

        # Time loop
        for i in range(len(time_obs)-1):

            # Prediction step (sample state transition)
            start_time = time.time()

            if is_parallel:
                if i == 0:
                    pool = mp.Pool(processes=mp.cpu_count(), initializer=self._set_local_propagator, initargs=())
                ret = pool.map(self._propagate_single, 
                    [(samples[i, j, :], time_obs[i], time_obs[i+1], GMe) for j in range(self.N)]
                )
                results = [r for r in ret]
                if i == len(time_obs)-2:
                    pool.close()
                # Postprocessing
                for j in range(self.N):
                    if results[j][0] == -1:
                        samples[i+1, j, :] = np.nan
                        weight[j] = 0.0
                        weight /= sum(weight)
                    else:
                        samples[i+1, j, :] = results[j][1]
            else:
                for j in range(self.N):  # loop over each particle
                    r, v = pk.par2ic(samples[i, j, 0:6], GMe)
                    prop = prop\
                        .state(np.concatenate((r, v)))\
                        .time(time_obs[i])\
                        .cram(samples[i, j, 6] * 1e-6)
                    ephem_j = prop([time_obs[i+1]])
                    # Handle infeasible solution 
                    if ephem_j.get_keplerian().size == 0:  # which means the particle is infeasible
                        samples[i+1, j, :] = np.nan
                        weight[j] = 0.0
                        weight /= sum(weight)
                    else:
                        samples[i+1, j, 0:6] = ephem_j.get_keplerian()[0]
                        samples[i+1, j, 6] = samples[i, j, 6]

            # Add process noise
            samples[i+1, :, :] += self.rng.multivariate_normal(np.zeros(7), self.Q, (self.N, ))
            end_time = time.time()
            print("Propagation of the particles is done. Wall-clock time: {} sec".format(end_time - start_time))

            # Weight update
            if i >= self.switch_R_after:
                self.R = self.R1
            likelihood = self._likelihood(obs_his[i, :], samples[i+1, :, :])
            weight = weight * likelihood
            weight[np.isnan(weight)] = 0.0
            weight /= sum(weight)

            # Resampling
            weight, states_out = self._resample_move(samples[i+1, :, :], weight, obs_his[i, :])
            samples[i+1, :, :] = states_out
            print("Cram at iter {0}: {1}".format(i, np.average(samples[i+1, :, 6]) * 1e-6))
            print("Maximum likelihood: {}".format(np.max(self._likelihood(obs_his[i, :], samples[i+1, :, :]))))
        
            # Compute state mean & covariance (diagonal)
            mean, cov_diag = self._calc_mean_cov(samples[i+1, :, :])
            mean_out[i+1, :] = mean
            cov_diag_out[i+1, :] = cov_diag

        return samples, mean_out, cov_diag_out

    def _set_local_propagator(self):
        global prop 
        prop = Propagator()

    def _propagate_single(self, params):
        global prop
        state, t0, tend, GMe = params
        r, v = pk.par2ic(state[0:6], GMe)
        prop = prop\
            .state(np.concatenate((r, v)))\
            .time(t0)\
            .cram(state[6] * 1e-6)
        ephem_j = prop([tend])
        # Handle infeasible solution 
        if ephem_j.get_keplerian().size == 0:  # which means the particle is infeasible
            flag = -1
            state_next = np.zeros(7)
        else:
            flag = 0
            state_next = np.append(ephem_j.get_keplerian()[0], state[6])
        return (flag, state_next)

    def _calc_mean_cov(self, state:np.ndarray):
        """
        state: np.ndarray whose size is (n_particles, 7)
        """
        diff_wrt0 = state - state[0, :]
        normalize_diff(diff_wrt0[:, 2])
        normalize_diff(diff_wrt0[:, 3])
        normalize_diff(diff_wrt0[:, 4])
        normalize_diff(diff_wrt0[:, 5])
        mean = state[0, :] + np.average(diff_wrt0, axis=0)
        diff_wrtavg = state - mean
        normalize_diff(diff_wrtavg[:, 2])
        normalize_diff(diff_wrtavg[:, 3])
        normalize_diff(diff_wrtavg[:, 4])
        normalize_diff(diff_wrtavg[:, 5])
        covariance = np.cov(diff_wrtavg, rowvar=False)
        return mean, np.diag(covariance)

    def _likelihood(self, obs:np.ndarray, state:np.ndarray):
        """
        obs: np.ndarray whose size is (n_obs, )
        state: np.ndarray whose size is (n_samples, 7)
        Returns likelihood as np.ndarray whose size is (n_samples, )
        """
        diff = obs - state[:, :6]
        normalize_diff(diff[:, 2])
        normalize_diff(diff[:, 3])
        normalize_diff(diff[:, 4])
        normalize_diff(diff[:, 5])
        mvn = multivariate_normal(np.zeros(6), self.R, allow_singular=True)
        likelihood = mvn.pdf(diff)
        return likelihood


    def _resample_move(self, states:np.ndarray, weights:np.ndarray, obs:np.ndarray):
        """
        state: np.ndarray whose size is (n_samples, 7)
        weights: np.ndarray whose size is (n_samples, )
        obs: np.ndarray whose size is (n_obs, )
        """
        # First, perform regular resampling
        # Compute cdf
        cdf_w = np.cumsum(weights)
        cdf_w_2d = np.tile(cdf_w, (self.N, 1)).transpose()  # (cdf axis, N)
        # Draw bunch of random #'s from U(0, 1)
        u_rands = self.rng.uniform(1e-10, 1.0, self.N)
        u_rands_2d = np.tile(u_rands, (self.N, 1))
        # Compare random numbers against CDF
        is_lower = u_rands_2d < cdf_w_2d
        is_lower_sum = np.cumsum(is_lower, axis=0)
        # Find a location where u < cdf occurs first
        ind, _ = np.nonzero(is_lower_sum == 1)
        states_out = states[ind, :]
        weights_out = 1/self.N * np.ones(self.N)

        # Perform MCMC move
        c = 0.5
        m_max = 10
        cram_avg = np.average(states_out[:, 6])
        cov = np.diag(np.array([50.0, 0.001, 0.002, 0.05, 0.05, 0.05, 0.02 * cram_avg])**2)
        for m in range(m_max):
            # Compute sample covariance
            # cov = np.cov(states_out, rowvar=False)
            # Add noise to the states
            e_m = self.rng.multivariate_normal(np.zeros(7), cov, (self.N, ))
            x_tilde = states_out + c * e_m
            # Compute acceptance ratio
            ratio = self._likelihood(obs, x_tilde) / self._likelihood(obs, states_out)
            # Acceptance/Rejection
            u_rands = self.rng.uniform(0.0, 1.0, self.N)
            is_accept = u_rands < ratio
            states_out[is_accept, :] = x_tilde[is_accept, :]
            acceptance_ratio = sum(is_accept)/len(is_accept)

        
        return weights_out, states_out

        


        
