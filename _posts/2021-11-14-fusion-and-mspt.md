---
layout: post
title: "Fusion reactors are a bit like molten salt power towers"
author: "Jacob Schwartz"
categories: journal
tags: [fusion,economics]
image: Crescent_Dunes_Solar_December_2014.jpg
use_math: false
---

*Image from [wikimedia.](https://en.wikipedia.org/wiki/Concentrated_solar_power#/media/File:Crescent_Dunes_Solar_December_2014.JPG)*
Fusion power plants don't yet exist, but it would be good to know what they *would* look like, especially in terms of cost and operational characteristics. I recently realized that molten salt [solar power towers](https://en.wikipedia.org/wiki/Solar_power_tower) (MSPTs) and more generally, [Concentrating Solar Plants](https://www.energy.gov/eere/solar/concentrating-solar-thermal-power) (CSPs), which currently do exist, have a lot in common with future fusion reactors. 

Molten salt power towers are a type of solar energy power plant which use a large array of mirrors to focus sunlight on a tank of hot liquid salt. The salt gets heated up to several hundred degrees Celsius. 
The hot salt can then be stored in a tank for up to several hours.
When the plant wants to generate electricity, the salt in the tank gets dispatched to boil water and make hot steam which then turns a turbine.
The ability to store the hot salt allows the plant to continue generating electricity after dark, which standard solar panels (photovoltaics) cannot. This makes them more flexible.

Here are six reasons why MSPTs are like future fusion reactors.

## High capital cost, low variable cost
In both plants, the heat generation system has a high upfront capital cost and low variable cost. For MSPTs the sunlight is free of course, and for fusion the cost of additional generation is very low.
It mostly comes down to wear and tear on the vacuum vessel and surrounding components, not fuel.
This is also similar to the situation in nuclear fission plants.

## Intermittent heat generation with fast transients
Large power plants like coal, natural gas, or nuclear can be limited in the speed that they ramp up and down the power.
For coal plants this may be limited by thermal expansion in the firebox, for nuclear power this can be limited by Xenon buildup effects, and for any thermal plants there are ramp rate limitations for the turbine and heat exchangers.

MSPTs and fusion reactors on the other hand may be required to deal with fast changes in initial heat generation.
Pulsed tokamaks would need to handle the thermal power from 0% to 100% and later from 100% to 0% in five minutes or less.
Even for steady state tokamaks the startup process is quick, it's not obvious to me that you can very slowly ramp up the power output.

MSPTs need to handle clouds! If a cloud blows over the site the power directed to the reciever will suddenly drop.
(Here it probably helps that the mirror field is a kilometer in diameter, so at least the change happens over a few minutes rather than instantly.) 

## Molten salt thermal energy storage
A bit of thermal buffer in the form of a large tank may be desirable for MSPTs in order to protect the steam generator from power transients when a cloud rolls over. (It must not be strictly required as some installations do not have dedicated storage.) 
Since the molten salts used are fairly dense with thermal energy and remain at low pressure, it's feasible to build a large tank to store them and extract the heat later in the day. This is already useful for MSPs as their power generation can be decoupled from that of solar panels.

Pulsed fusion reactors also may need a relatively small molten salt 'buffer' tank to protect the heat exchangers and generator from the pulsed behavior of the heat source. This is being studied for [EU-DEMO](https://linkinghub.elsevier.com/retrieve/pii/S0920379621002805).
Since a fusion core has a high capital cost and a low variable cost, we probably want to have it running as much as possible once it's built. This could be accomplished by extending the tank size in order to store up heat during the day and release it to generate electricity at night when solar PV is unavailable.

## Unit size

The largest solar power towers generate roughly 100-150MWe. This is similar to the sizes we might expect for early fusion reactors, and is an attractive size for power companies.

## Potential for advanced energy conversion

Higher temperature leads to the possibility of more efficient power conversion. While current MSPTs use nitrate salt mixture with a maximum temperature of about 565°C, researchers working on 'Gen3 Concentrating Solar Power' technologies are looking at using materials that can run hotter, like 700°C. Then it would be attractive to utilize power cycles with supercitical CO₂ or He, which may have higher effienciency, smaller physical size, and quicker ramping capabilities.

Fusion reactors are limited in their working fluid temperature only by our engineering abilities, not by the reaction itself. Engineers are also examining and designing novel power cycles to work with fusion.
Since fusion reactors often require significant recirculating electric power, there's an especially large incentive to have a highly efficient power cycle.

## Extreme heat fluxes

(This is more of a physics similarity than an engineering similarity perhaps.)
Fusion reactors and MSPTs both need to push heat through a limited surface area. 
Unlike in *fission* reactors, where the entire core is essentially one large heat exchanger, (tokamak) fusion reactors need to push the heat they generate through the vacuum vessel wall before it reaches the working fluid. Even machines that directly employ a liquid metal wall (Z-pinches or magneto-inertial fusion) may need to be careful not to flash-evaporate the liquid surface.

MSPTs focus sunlight on the 'reciever' atop the tower, and the light reaches heat fluxes of several megawatts per square meter. Fusion reactor walls see similar and larger heat fluxes.
(Ultimately the technology here is not directly transferrable as fusion devices operate in high vacuum and radiation environments. MSPTs operate in air with an oxidizing environment.)

# Key takeaway

As far as I can tell there's not a *huge* amount of current research on developing sCO₂ or He Brayton power cycles, nor on technology to handle high-temperature molten salts. (Though see these recent [awards.](https://www.solarpaces.org/these-25-advanced-csp-cst-technologies-to-share-33-million-us-doe-funding/))
Perhaps Gen3 Concentrating Solar Power researchers and fusion engineers should collaborate somehow to develop these advanced power cycle technologies, which would be mutually beneficial.
