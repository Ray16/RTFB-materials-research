"""Reservations, host capacity, ownership and checkpoint transport.

GPU work enters through ``gpu_reserve``. This package does not choose a device
by setting ``CUDA_VISIBLE_DEVICES`` itself.
"""
