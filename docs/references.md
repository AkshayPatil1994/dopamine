---
title: References
nav_order: 10
---

<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# References

Full bibliography of the works cited in these pages.

<ol>
{% for item in site.data.references %}
  {% assign r = item[1] %}
  <li id="{{ item[0] }}">{{ r.authors }}, "{{ r.title }}," <em>{{ r.venue }}</em>, {{ r.year }}.
  {% if r.doi %}<a href="https://doi.org/{{ r.doi }}">doi:{{ r.doi }}</a>{% endif %}</li>
{% endfor %}
</ol>
