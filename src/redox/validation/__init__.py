"""Compatibility paths for the certified occupied-sampling snapshot producer.

``quality/sampling_snapshot.py`` is byte-identical to the producer named by the
frozen receipt audit. That file imports these historical module paths. New code
imports ``redox.quality`` directly.
"""
